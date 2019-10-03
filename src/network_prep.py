import datetime
import os
import sys
import time
import yaml
import pandas
import redis
from knpackage.toolbox import get_run_parameters, get_run_directory_and_file
import utils.log_util as logger
from utils.io_util import IOUtil

MGET_CHUNK = 10000
KEEP_THR = 0.6

def load_edge_file(file_path):
    """
    Loads edge file as a DataFrame object.

    Args:
        file_path: input file, which is uploaded from frontend

    Returns:
        input_df: user input as a DataFrame, which doesn't have any header or index
    """
    if not file_path or not file_path.strip() or not os.path.exists(file_path):
        logger.logging.append('ERROR: Input file path is not valid: {}. Please provide a valid '
                              'input path.'.format(file_path))
        return None
    try:
        # loads input data
        input_df = pandas.read_csv(file_path, sep='\t', header=None, index_col=None)
        if input_df.shape == (0, 0):
            logger.logging.append('ERROR: Input data {} is empty. Please provide a valid '
                                  'input data.'.format(file_path))
            return None
        if input_df.shape[1] < 3:
            logger.logging.append('ERROR: Not enough columns in input data {}. Please provide '
                                  'a valid input data.'.format(file_path))
            return None
        logger.logging.append(
            'INFO: Successfully loaded input data: {} with {} row(s) and {} '
            'column(s)'.format(file_path, input_df.shape[0], input_df.shape[1]))
        return input_df
    except Exception as err:
        logger.logging.append('ERROR: {}'.format(str(err)))
        return None

def check_column_data(dataframe, col_idx=2, chk_na=False, chk_real_num=False, chk_pos_num=False):
    """
    Customized checks for input data (contains NA value, all real numbers, all positive numbers)
    Args:
        dataframe: input DataFrame to be checked
        col_idx: column number (from 0) to be checked for proper values
        chk_na: check NA in DataFrame column
        chk_real_num: check only real number exists in DataFrame column
        chk_pos_num: check only positive number exists in DataFrame column
    Returns:
        dataframe: cleaned DataFrame
    """
    # checks if dataframe contains NA value
    if chk_na is True:
        if dataframe[[col_idx]].isnull().values.any():
            logger.logging.append("ERROR: Found NaN in value column " + str(col_idx+1))
            return None
    # checks real number negative to positive infinite
    if chk_real_num is True:
        if False in dataframe[[col_idx]].applymap(lambda x: isinstance(x, (int, float))).values:
            logger.logging.append("ERROR: Found non-numeric value in column " + str(col_idx+1))
            return None
    # checks if dataframe contains only non-negative number
    if chk_pos_num is True:
        if False in dataframe[[col_idx]].applymap(lambda x: x >= 0).values:
            logger.logging.append("ERROR: Found negative value in column " + str(col_idx+1))
            return None
    return dataframe

def get_database(redis_host, redis_port, redis_pass):
    """Returns a Redis database connection.
    This returns a Redis database connection access to its functions if the
    module is imported.
    Args:
        args (Namespace): args as populated namespace or 'None' for defaults
    Returns:
        StrictRedis: a redis connection object
    """
    return redis.StrictRedis(host=redis_host, port=redis_port,
                             password=redis_pass)

def get_node_info(rdb, fk_array, ntype, hint, taxid):
    """Uses the redis database to convert a node alias to KN internal id

    Figures out the type of node for each id in fk_array and then returns
    all of the metadata associated or unmapped-*

    Args:
        rdb (redis object): redis connection to the mapping db
        fk_array (list): the array of foreign gene identifers to be translated
        ntype (str): 'Gene' or 'Property' or None
        hint (str): a hint for conversion
        taxid (str): the species taxid, None if unknown

    Returns:
        list: list of lists containing 5 col info for each mapped gene
    """
    hint = None if hint == '' or hint is None else hint.upper()
    taxid = None if taxid == '' or taxid is None else str(taxid)
    if ntype == '':
        ntype = None

    if ntype is None:
        res_arr = rdb.mget(['::'.join(['stable', str(fk), 'type']) for fk in fk_array])
        fk_prop = [fk for fk, res in zip(fk_array, res_arr) if res is not None
                   and res.decode() == 'Property']
        fk_gene = [fk for fk, res in zip(fk_array, res_arr) if res is not None
                   and res.decode() == 'Gene']
        if fk_prop and fk_gene:
            raise ValueError("Mixture of property and gene nodes.")
        ntype = 'Property' if fk_prop else 'Gene'

    if ntype == "Gene":
        stable_array = conv_gene(rdb, fk_array, hint, taxid)
    elif ntype == "Property":
        stable_array = fk_array
    else:
        raise ValueError("Invalid ntype")

    return list(zip(fk_array, *node_desc(rdb, stable_array)))

def conv_gene(rdb, fk_array, hint, taxid):
    """Uses the redis database to convert a gene to ensembl stable id

    This checks first if there is a unique name for the provided foreign key.
    If not it uses the hint and taxid to try and filter the foreign key
    possiblities to find a matching stable id.

    Args:
        rdb (redis object): redis connection to the mapping db
        fk_array (list): the foreign gene identifers to be translated
        hint (str): a hint for conversion
        taxid (str): the species taxid, 'unknown' if unknown

    Returns:
        str: result of searching for gene in redis DB
    """
    hint = None if hint == '' or hint is None else hint.upper()
    taxid = None if taxid == '' or taxid is None else str(taxid)

    #use ensembl internal uniprot mappings
    if hint == 'UNIPROT' or hint == 'UNIPROTKB':
        hint = 'UNIPROT_GN'

    ret_stable = ['unmapped-none'] * len(fk_array)

    def replace_none(ret_st, pattern):
        """Search redis for genes that still are unmapped
        """
        curr_none = [i for i in range(len(fk_array)) if ret_st[i] == 'unmapped-none']
        while curr_none:
            temp_curr_none = curr_none[:MGET_CHUNK]
            curr_none = curr_none[MGET_CHUNK:]
            vals_array = rdb.mget([pattern.format(str(fk_array[i]).upper(), taxid, hint)
                                   for i in temp_curr_none])
            for i, val in zip(temp_curr_none, vals_array):
                if val is None:
                    continue
                ret_st[i] = val.decode()

    if hint is not None and taxid is not None:
        replace_none(ret_stable, 'triplet::{0}::{1}::{2}')
    if taxid is not None:
        replace_none(ret_stable, 'taxon::{0}::{1}')
    if hint is not None:
        replace_none(ret_stable, 'hint::{0}::{2}')
    if taxid is None:
        replace_none(ret_stable, 'unique::{0}')
    return ret_stable

def node_desc(rdb, stable_array):
    """Uses the redis database to find metadata about node given its stable id

    Return all metadata for each element of stable_array

    Args:
        rdb (redis object): redis connection to the mapping db
        stable_array (str): the array of stable identifers to be searched

    Returns:
        list: list of lists containing 4 col info for each mapped node
    """
    ret_type = ["None"] * len(stable_array)
    ret_alias = ["None"] * len(stable_array)
    ret_desc = ["None"] * len(stable_array)
    ret_biotype = ["None"] * len(stable_array)
    st_map_idxs = [idx for idx, st in enumerate(stable_array) if not st.startswith('unmapped')]

    while st_map_idxs:
        temp_idxs = st_map_idxs[:MGET_CHUNK]
        st_map_idxs = st_map_idxs[MGET_CHUNK:]

        type_array = rdb.mget(['::'.join(['stable', stable_array[i], 'type'])
                               for i in temp_idxs])
        alias_array = rdb.mget(['::'.join(['stable', stable_array[i], 'alias'])
                                for i in temp_idxs])
        desc_array = rdb.mget(['::'.join(['stable', stable_array[i], 'desc'])
                               for i in temp_idxs])
        btype_array = rdb.mget(['::'.join(['stable', stable_array[i], 'biotype'])
                                for i in temp_idxs])

        for i, val in zip(temp_idxs, type_array):
            if val is None:
                continue
            ret_type[i] = val.decode()
        for i, val in zip(temp_idxs, alias_array):
            if val is None:
                continue
            ret_alias[i] = val.decode()
        for i, val in zip(temp_idxs, desc_array):
            if val is None:
                continue
            ret_desc[i] = val.decode()
        for i, val in zip(temp_idxs, btype_array):
            if val is None:
                continue
            ret_biotype[i] = val.decode()

    return stable_array, ret_type, ret_alias, ret_desc, ret_biotype


class Pipelines:
    def __init__(self, run_parameters):
        self.run_parameters = run_parameters
        if 'network_threshold' not in self.run_parameters.keys():
            self.run_parameters['network_threshold'] = KEEP_THR
        if 'make_symmetric' not in self.run_parameters.keys():
            self.run_parameters['network_threshold'] = True
        # load edge file
        self.raw_network_df = load_edge_file(
            self.run_parameters['raw_edgefile_full_path']) \
            if 'raw_edgefile_full_path' in self.run_parameters.keys() else None
    def network_prepper_pipeline(self):
        """
        Runs data cleaning for network mapping pipeline.
        Returns:
            validation_flag: Boolean type value indicating if input data is valid or not.
            message: A message indicates the status of current check.
        """
        if self.raw_network_df is None:
            return False, logger.logging

        # Checks only non-negative real number appears in column three
        raw_network_chked = check_column_data(self.raw_network_df, 2,
                                              chk_na=True, chk_real_num=True,
                                              chk_pos_num=True)
        if raw_network_chked is None:
            return False, logger.logging

        logger.logging.append('INFO: Checked raw edgefile has {} row(s), {} column(s).'.format(
            raw_network_chked.shape[0], raw_network_chked.shape[1]))

        # maps the first two node columns
        try:
            rdb = get_database(self.run_parameters['redis_credential']['host'],
                               self.run_parameters['redis_credential']['port'],
                               self.run_parameters['redis_credential']['password'])
            sources_mapped = get_node_info(rdb, raw_network_chked.iloc[:, 0].tolist(),
                                           "Gene", self.run_parameters['source_hint'],
                                           self.run_parameters['taxonid'])
            targets_mapped = get_node_info(rdb, raw_network_chked.iloc[:, 1].tolist(),
                                           "Gene", self.run_parameters['source_hint'],
                                           self.run_parameters['taxonid'])
        except Exception as err:
            logger.logging.append('ERROR: Problem with Redis mapping: {} : {}'.format(
                self.run_parameters['redis_credential']['host'], str(err)))
            return False, logger.logging

        # prints out the full edge file with mapping and the full node mapping
        try:
            mapped_sources_df = pandas.DataFrame(sources_mapped)
            mapped_targets_df = pandas.DataFrame(targets_mapped)

            mapped_edge_df = pandas.concat([mapped_sources_df[[1, 3]],
                                            mapped_targets_df[[1, 3]],
                                            raw_network_chked], axis=1)
            IOUtil.write_to_file(mapped_edge_df, self.run_parameters['raw_edgefile_full_path'],
                                 self.run_parameters['results_directory'], '.full_mapped_edges.tsv',
                                 use_index=False, use_header=False)
            logger.logging.append('INFO: Full mapped edgefile has {} row(s), {} column(s).'.format(
                mapped_edge_df.shape[0], mapped_edge_df.shape[1]))

            mapped_nodes_df = pandas.concat([mapped_sources_df, mapped_targets_df])
            mapped_nodes_df.sort_values(mapped_nodes_df.columns[0], inplace=True)
            mapped_nodes_df.drop_duplicates(keep='first', inplace=True)
            IOUtil.write_to_file(mapped_nodes_df, self.run_parameters['raw_edgefile_full_path'],
                                 self.run_parameters['results_directory'], '.full_mapped_nodes.tsv',
                                 use_index=False, use_header=False)
            logger.logging.append('INFO: Full mapped nodefile has {} row(s), {} column(s).'.format(
                mapped_nodes_df.shape[0], mapped_nodes_df.shape[1]))
        except Exception as err:
            logger.logging.append('ERROR: Problem with printing full mapped results: {}'.format(
                str(err)))
            return False, logger.logging

        # print out clean/deduplicated edge and node files
        # drop Nones
        clean_nodes_df = mapped_nodes_df[mapped_nodes_df.iloc[:, 2] != "None"]
        logger.logging.append('INFO: {} mapped nodes out of {} original.'.format(
            clean_nodes_df.shape[0], mapped_nodes_df.shape[0]))
        if ((clean_nodes_df.shape[0] / mapped_nodes_df.shape[0]) <
                self.run_parameters['network_threshold']):
            logger.logging.append('ERROR: Less than {}% of nodes are mapped (only {} '
                                  'out of {}.)'.format(
                                      self.run_parameters['network_threshold']*100,
                                      clean_nodes_df.shape[0],
                                      mapped_nodes_df.shape[0]))
            return False, logger.logging

        clean_edge_df = mapped_edge_df.copy()
        clean_edge_df = clean_edge_df[clean_edge_df.iloc[:, 1] != "None"]
        clean_edge_df = clean_edge_df[clean_edge_df.iloc[:, 3] != "None"]
        clean_edge_df = clean_edge_df.iloc[:, [0, 2, 6]].reset_index(drop=True)
        clean_edge_df.columns = ['source', 'target', 'score']
        logger.logging.append('INFO: {} mapped edges out of {} original.'.format(
            clean_edge_df.shape[0], mapped_edge_df.shape[0]))
        if ((clean_edge_df.shape[0] / mapped_edge_df.shape[0]) <
                self.run_parameters['network_threshold']):
            logger.logging.append('ERROR: Less than {} of edges have both nodes '
                                  'mapped (only {} out of {}).'.format(
                                      self.run_parameters['network_threshold']*100,
                                      clean_edge_df.shape[0],
                                      mapped_edge_df.shape[0]))
            return False, logger.logging

        # make symmetric
        if self.run_parameters['make_symmetric']:
            df2 = clean_edge_df.copy().reset_index(drop=True)
            df2.columns = ['target', 'source', 'score']
            clean_edge_df = pandas.concat([clean_edge_df, df2])[['source', 'target', 'score']]
            logger.logging.append('INFO: mapped edges were made symmetric.')
        # sort by name and score
        clean_edge_df.sort_values(by=['source', 'target', 'score'],
                                  inplace=True, ascending=False)
        # deduplicate
        clean_edge_df.drop_duplicates(subset=['source', 'target'], keep='first', inplace=True)
        logger.logging.append('INFO: {} mapped edges after duplicates removed.'.format(
            clean_edge_df.shape[0]))

        # make unweighted
        # normalize edge weights
        # make upper triangle

        # check if fewer edges than a tree
        if self.run_parameters['make_symmetric']:
            if clean_edge_df.shape[0] < 2*(clean_nodes_df.shape[0] - 1):
                logger.logging.append('ERROR: Clean graph is too sparce with only {} symmetric '
                                      'edges for {} clean nodes.'.format(clean_edge_df.shape[0],
                                                                         clean_nodes_df.shape[0]))
                return False, logger.logging
        else:
            if clean_edge_df.shape[0] < (clean_nodes_df.shape[0] - 1):
                logger.logging.append('ERROR: Clean graph is too sparce with only {} clean edges '
                                      'for {} clean nodes.'.format(clean_edge_df.shape[0],
                                                                   clean_nodes_df.shape[0]))
                return False, logger.logging

        # write clean files
        IOUtil.write_to_file(clean_nodes_df, self.run_parameters['raw_edgefile_full_path'],
                             self.run_parameters['results_directory'], '.clean.node_map',
                             use_index=False, use_header=False)
        logger.logging.append('INFO: Clean node_map file has {} row(s), {} column(s).'.format(
            clean_nodes_df.shape[0], clean_nodes_df.shape[1]))

        IOUtil.write_to_file(clean_edge_df, self.run_parameters['raw_edgefile_full_path'],
                             self.run_parameters['results_directory'], '.clean.edge',
                             use_index=False, use_header=False)
        logger.logging.append('INFO: Clean edge file has {} row(s), {} column(s).'.format(
            clean_edge_df.shape[0], clean_edge_df.shape[1]))

        output_file_basename = os.path.splitext(os.path.basename(
            os.path.normpath(self.run_parameters['raw_edgefile_full_path'])))[0]
        total_edges = clean_nodes_df.shape[0] * clean_nodes_df.shape[0]
        dict_file = {
            'data': {
                'density': clean_edge_df.shape[0] / total_edges,
                'num_connected_components': -1,
                'num_edges': clean_edge_df.shape[0],
                'num_gene_nodes': clean_nodes_df.shape[0],
                'num_nodes': clean_nodes_df.shape[0],
                'num_prop_nodes': 0
            },
            'raw_stats': {
                'num_orig_edges': mapped_edge_df.shape[0],
                'num_orig_nodes': mapped_nodes_df.shape[0]
            },
            "id": output_file_basename,
            "datasets": {
                'orig_file': self.run_parameters['raw_edgefile_full_path'],
                'conv_date': datetime.datetime.utcfromtimestamp(int(time.time()))
            },
            "species": {
                "taxon_identifier": self.run_parameters['taxonid']
            },
            "edge_type": {
                "n1_type": "Gene",
                "n2_type": "Gene",
                "type_desc": "custom uploaded network",
                "score_desc": "custom score uploaded",
                "score_best": clean_edge_df['score'].max().item(),
                "score_worst": clean_edge_df['score'].min().item(),
                "bidirectional": 1
            },
            'run_params': self.run_parameters
        }
        metayaml = (self.run_parameters['results_directory'] + '/' +
                    output_file_basename + ".metadata")
        with open(metayaml, 'w') as file:
            yaml.dump(dict_file, file, default_flow_style=False)

        return True, logger.logging



def run_pipelines(run_parameters, method):
    validation_flag, message = getattr(Pipelines(run_parameters), method)()
    logger.generate_logging(validation_flag, message,
                            run_parameters["results_directory"] + "/log_" + method + ".yml")


def network_prepper():
    try:
        logger.init()
        run_directory, run_file = get_run_directory_and_file(sys.argv)
        run_parameters = get_run_parameters(run_directory, run_file)
        run_pipelines(run_parameters, "network_prepper_pipeline")
    except Exception as err:
        logger.logging.append("ERROR: {}".format(str(err)))
        # try to write the log
        logger.generate_logging(False, logger.logging,
                                run_parameters["results_directory"] + "/log_" +
                                "network_prepper_pipeline" + ".yml")
        raise RuntimeError(str(err))


if __name__ == "__main__":
    network_prepper()
