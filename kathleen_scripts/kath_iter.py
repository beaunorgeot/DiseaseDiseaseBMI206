# iterator for separation.py

# import modules

import itertools
import os, csv, glob

listing = glob.glob('disease_csvs/*.txt')

for c, v in itertools.combinations(listing,2):
    import networkx as nx
    import numpy as np
    import optparse
    import sys, os

    '''
    # =============================================================================

               S T A R T   D E F I N I T I O N S 

    # =============================================================================
    '''
    def read_network(network_file):
        """
        Reads a network from an external file.

        * The edgelist must be provided as a tab-separated table. The
        first two columns of the table will be interpreted as an
        interaction gene1 <==> gene2

        * Lines that start with '#' will be ignored
        """

        G = nx.Graph()
        for line in open(network_file,'r'):
            # lines starting with '#' will be ignored
            if line[0]=='#':
                continue
            # The first two columns in the line will be interpreted as an
            # interaction gene1 <=> gene2
            line_data   = line.strip().split('\t')
            node1 = line_data[0]
            node2 = line_data[1]
            G.add_edge(node1,node2)

        print "\n> done loading network:"
        print "> network contains %s nodes and %s links" %(G.number_of_nodes(),
                                                           G.number_of_edges())
        
        return G


    # =============================================================================
    def read_gene_list(gene_file):
        """
        Reads a list genes from an external file.

        * The genes must be provided as a table. If the table has more
        than one column, they must be tab-separated. The first column will
        be used only.

        * Lines that start with '#' will be ignored
        """

        genes_set = set()
        for line in open(gene_file,'r'):
            # lines starting with '#' will be ignored
            if line[0]=='#':
                continue
            # the first column in the line will be interpreted as a seed
            # gene:
            line_data = line.strip().split('\t')
            gene      = line_data[0]
            genes_set.add(gene)

        print "\n> done reading genes:"
        print "> %s genes found in %s" %(len(genes_set),gene_file)

        return genes_set


    # =============================================================================
    def remove_self_links(G):

        sl = G.selfloop_edges()
        G.remove_edges_from(sl)



    # =============================================================================
    def get_pathlengths_for_single_set(G,given_gene_set):
        
        """
        calculate the shortest paths of a given set of genes in a
        given network. The results are stored in a dictionary of
        dictionaries:
        all_path_lenghts[gene1][gene2] = l
        with gene1 < gene2, so each pair is stored only once!

        PARAMETERS:
        -----------
            - G: network
            - gene_set: gene set for which paths should be computed

        RETURNS:
        --------
            - all_path_lenghts[gene1][gene2] = l for all pairs of genes
              with gene1 < gene2

        """ 

        # remove all nodes that are not in the network
        all_genes_in_network = set(G.nodes())
        gene_set = given_gene_set & all_genes_in_network

        all_path_lenghts = {}
        
        # calculate the distance of all possible pairs
        for gene1 in gene_set:
            if not all_path_lenghts.has_key(gene1):
                all_path_lenghts[gene1] = {}
            for gene2 in gene_set:
                if gene1 < gene2:
                    try:
                        l = nx.shortest_path_length(G, source=gene1, target=gene2)
                        all_path_lenghts[gene1][gene2] = l
                    except:
                        continue

        return all_path_lenghts



    # =============================================================================
    def get_pathlengths_for_two_sets(G,given_gene_set1,given_gene_set2):
        
        """
        calculate the shortest paths between two given set of genes in a
        given network. The results are stored in a dictionary of
        dictionaries: all_path_lenghts[gene1][gene2] = l with gene1 <
        gene2, so each pair is stored only once!

        PARAMETERS:
        -----------
            - G: network
            - gene_set1/2: gene sets for which paths should be computed

        RETURNS:
        --------
            - all_path_lenghts[gene1][gene2] = l for all pairs of genes
              with gene1 < gene2

        """ 

        # remove all nodes that are not in the network
        all_genes_in_network = set(G.nodes())
        gene_set1 = given_gene_set1 & all_genes_in_network
        gene_set2 = given_gene_set2 & all_genes_in_network

        all_path_lenghts = {}
        
        # calculate the distance of all possible pairs
        for gene1 in gene_set1:
            if not all_path_lenghts.has_key(gene1):
                all_path_lenghts[gene1] = {}
            for gene2 in gene_set2:
                if gene1 != gene2:
                    try:
                        l = nx.shortest_path_length(G, source=gene1, target=gene2)
                        if gene1 < gene2:
                            all_path_lenghts[gene1][gene2] = l
                        else:
                            if not all_path_lenghts.has_key(gene2):
                                all_path_lenghts[gene2] = {}
                            all_path_lenghts[gene2][gene1] = l
                    except:
                        continue

        return all_path_lenghts


    # =============================================================================
    def calc_single_set_distance(G,given_gene_set):

        """
        Calculates the mean shortest distance for a set of genes on a
        given network    
        

        PARAMETERS:
        -----------
            - G: network
            - gene_set: gene set for which distance will be computed 

        RETURNS:
        --------
             - mean shortest distance 

        """


        # remove all nodes that are not in the network, just to be safe
        all_genes_in_network = set(G.nodes())
        gene_set = given_gene_set & all_genes_in_network

        # get the network distances for all gene pairs:
        all_path_lenghts = get_pathlengths_for_single_set(G,gene_set)

        all_distances = []

        # going through all gene pairs
        for geneA in gene_set:

            all_distances_A = []
            for geneB in gene_set:

                # I have to check which gene is 'smaller' in order to know
                # where to look up the distance of that pair in the
                # all_path_lengths dict
                if geneA < geneB:
                    if all_path_lenghts[geneA].has_key(geneB):
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                else:
                    if all_path_lenghts[geneB].has_key(geneA):
                        all_distances_A.append(all_path_lenghts[geneB][geneA])

            if len(all_distances_A) > 0:
                l_min = min(all_distances_A)
                all_distances.append(l_min)

        # calculate mean shortest distance
        mean_shortest_distance = np.mean(all_distances)

        return mean_shortest_distance


    # =============================================================================
    def calc_set_pair_distances(G,given_gene_set1,given_gene_set2):

        """
        Calculates the mean shortest distance between two sets of genes on
        a given network
        
        PARAMETERS:
        -----------
            - G: network
            - gene_set1/2: gene sets for which distance will be computed 

        RETURNS:
        --------
             - mean shortest distance 

        """

        # remove all nodes that are not in the network
        all_genes_in_network = set(G.nodes())
        gene_set1 = given_gene_set1 & all_genes_in_network
        gene_set2 = given_gene_set2 & all_genes_in_network

        # get the network distances for all gene pairs:
        all_path_lenghts = get_pathlengths_for_two_sets(G,gene_set1,gene_set2)

        all_distances = []

        # going through all pairs starting from set 1 
        for geneA in gene_set1:

            all_distances_A = []
            for geneB in gene_set2:

                # the genes are the same, so their distance is 0
                if geneA == geneB:
                    all_distances_A.append(0)
                    
                # I have to check which gene is 'smaller' in order to know
                # where to look up the distance of that pair in the
                # all_path_lengths dict
                else:
                    if geneA < geneB:
                        try:
                            all_distances_A.append(all_path_lenghts[geneA][geneB])
                        except:
                            pass

                    else:
                        try:
                            all_distances_A.append(all_path_lenghts[geneB][geneA])
                        except:
                            pass


            if len(all_distances_A) > 0:
                l_min = min(all_distances_A)
                all_distances.append(l_min)

        # going through all pairs starting from disease B
        for geneA in gene_set2:

            all_distances_A = []
            for geneB in gene_set1:

                # the genes are the same, so their distance is 0
                if geneA == geneB:
                    all_distances_A.append(0)

                # I have to check which gene is 'smaller' in order to know
                # where to look up the distance of that pair in the
                # all_path_lengths dict
                else:
                    if geneA < geneB:
                        try:
                            all_distances_A.append(all_path_lenghts[geneA][geneB])
                        except:
                            pass
                            
                    else:
                        try:
                            all_distances_A.append(all_path_lenghts[geneB][geneA])
                        except:
                            pass

            if len(all_distances_A) > 0:
                l_min = min(all_distances_A)
                all_distances.append(l_min)


        # calculate mean shortest distance
        mean_shortest_distance = np.mean(all_distances)

        return mean_shortest_distance



    """
    # =============================================================================

               E N D    O F    D E F I N I T I O N S 

    # =============================================================================
    """


    if __name__ == '__main__':

        # "Hey Ho, Let's go!" -- The Ramones (1976)

        # --------------------------------------------------------
        # 
        # PARSING THE COMMAND LINE
        # 
        # --------------------------------------------------------

        # parser = optparse.OptionParser()

        # parser.add_option('-u', '--usage',
        #                   help    ='print more info on how to use this script',
        #                   action="callback", callback=print_usage)

        # parser.add_option('-n',
        #                   help    ='file containing the network edgelist [interactome.tsv]',
        #                   dest    ='network_file',
        #                   default ='interactome.tsv',
        #                   type    = "string")

        # parser.add_option('--g1',
        #                   help    ='file containing gene set 1',
        #                   dest    ='gene_file_1',
        #                   default ='none',
        #                   type    = "string")

        # parser.add_option('--g2',
        #                   help    ='file containing gene set 2',
        #                   dest    ='gene_file_2',
        #                   default ='none',
        #                   type    = "string")

        # parser.add_option('-o',
        #                   help    ='file for results [separation_results.txt]',
        #                   dest    ='results_file',
        #                   default ='separation_results.txt',
        #                   type    = "string")


        # (opts, args) = parser.parse_args()

        network_file = 'interactome.tsv'
        gene_file_1  = c
        gene_file_2  = v

        # checking for input:
        if gene_file_1 == 'none' or gene_file_2 == 'none':
            error_message = """
            ERROR: you must specify two files with gene sets, for example:
            ./separation.py --g1 MS.txt --g2 PD.txt

            For more information, type
            ./separation.py --usage
            
            """
            print error_message
            sys.exit(0)

        if network_file == 'interactome.tsv':
            print '> default network from "interactome.tsv" will be used'


        # --------------------------------------------------------
        #
        # LOADING NETWORK and DISEASE GENES
        #
        # --------------------------------------------------------

        # read network
        G  = read_network(network_file)
        # get all genes ad remove self links
        all_genes_in_network = set(G.nodes())
        remove_self_links(G)

        # read gene set 1
        genes_A_full = read_gene_list(gene_file_1)
        # removing genes that are not in the network:
        genes_A = genes_A_full & all_genes_in_network
        if len(genes_A_full) != len(genes_A):
            print "> ignoring %s genes that are not in the network" %(
                len(genes_A_full - all_genes_in_network))
            print "> remaining number of genes: %s" %(len(genes_A))


        # read gene set 1
        genes_B_full = read_gene_list(gene_file_2)
        # removing genes that are not in the network:
        genes_B = genes_B_full & all_genes_in_network
        if len(genes_B_full) != len(genes_B):
            print "> ignoring %s genes that are not in the network" %(
                len(genes_B_full - all_genes_in_network))
            print "> remaining number of genes: %s" %(len(genes_B))


        # --------------------------------------------------------
        #
        # CALCULATE NETWORK QUANTITIES
        #
        # --------------------------------------------------------

        # distances WITHIN the two gene sets:
        d_A = calc_single_set_distance(G,genes_A)
        d_B = calc_single_set_distance(G,genes_B)

        # distances BETWEEN the two gene sets:
        d_AB = calc_set_pair_distances(G,genes_A,genes_B)

        # calculate separation
        s_AB = d_AB - (d_A + d_B)/2.

        with open('output_distances.tsv', 'ab') as f:
            row = str(c[13:]) +'\t'+ str(v[13:]) +'\t' + str(d_AB) +'\t'+ str(s_AB) + '\n'
            f.write((row))

    # for c in itertools.combinations(os.listdir('/disease_csvs/',2):
