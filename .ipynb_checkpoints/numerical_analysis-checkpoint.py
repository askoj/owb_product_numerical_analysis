# Import the required modules
import itertools
import math
import numpy as np
import decimal
decimal.getcontext().prec = 100
import matplotlib.pyplot as plt

'''
	This function produces the number of hyperedges within the 'Obeid-Wolfe-Bruza' (OWB) product.
'''
def prod_OWB_edges(num_parties, num_edges, num_outcomes):
	y = []
	# Compute the combinations of the various parties
	for i in range(num_parties-1, 0, -1):
		y.extend([x for x in itertools.combinations(list(range(num_parties)),i)])
	total_edges = 0
	# For each combination and its symmetric difference on the total parties, procure the edges
	# for its respective relation.
	for p in y:
		num_edges_signaller = (num_edges**(len(p)))*((num_outcomes-1)**(len(p)))
		num_edges_signalled = (num_edges**(num_parties-len(p)))-1
		total_add = num_edges_signaller*num_edges_signalled
		total_edges += total_add
	# Add on the hyperedges from the Cartesian product, and return the result
	return (total_edges+(num_edges**num_parties))

'''
	This function returns the hyperedges that constitute the theoretical minimum number of hyperedges,
	as related by Dr Elie Wolfe, and derived from his and Sainz' 'rref' method within 
	Sainz, A. B. and Wolfe, E. (2018). 
	Multipartite composition of contextuality scenarios. Foundations of Physics, 48(8):925–953.
'''
def prod_Elie_edges(num_parties, num_edges, num_outcomes):
	return (((num_outcomes*num_edges)**num_parties) - (((num_edges*(num_outcomes-1)) + 1)**num_parties) + 1)

'''
	This function produces the number of hyperedges within the minimal FR product
	Acin, A., Fritz, T., Leverrier, A., and Sainz, A. B. (2015). A combinatorial approach to nonlocality and contextuality. Communications in Mathematical Physics, 334(2):533–628.
'''
def prod_minFR_edges(num_parties, num_edges, num_outcomes):
	return ((num_parties*(num_edges**(num_parties-1)))*(num_edges**(num_outcomes**(num_parties-1)))) - ((num_parties-1)*(num_edges**num_parties))

'''
	This function produces the number of hyperedges within the common FR product
	Acin, A., Fritz, T., Leverrier, A., and Sainz, A. B. (2015). A combinatorial approach to nonlocality and contextuality. Communications in Mathematical Physics, 334(2):533–628.
'''
def prod_commonFR_edges(num_parties, num_edges, num_outcomes):
	parties = decimal.Decimal(num_parties)
	outcomes = decimal.Decimal(num_outcomes)
	edges = decimal.Decimal(num_edges)

	def AotimesB(edges_of_A, edges_of_B, outcomes_of_one_edge_of_A, outcomes_of_one_edge_of_B):
		#number of edges in A \otimes B
		edges_of_AB = (((edges_of_B**outcomes_of_one_edge_of_A))*edges_of_A)
		edges_of_BA = (((edges_of_A**outcomes_of_one_edge_of_B))*edges_of_B)
		#number of outcomes in any edge of A \otimes B
		outcomes_of_one_edge_of_AB = outcomes_of_one_edge_of_A*outcomes_of_one_edge_of_B
		return edges_of_AB+edges_of_BA, outcomes_of_one_edge_of_AB

	outcomes_in_composite_party = outcomes
	edges_in_composite_party = edges
	for i in range(int(parties)-1):
		edges_in_composite_party, outcomes_in_composite_party = AotimesB(edges, edges_in_composite_party, outcomes, outcomes_in_composite_party)
	#print(edges_in_composite_party+(edges**parties))
	if (num_parties > 8):
		return None
	return int((edges_in_composite_party)-(edges**parties))

'''
	This function produces the number of hyperedges within the maximal FR product
'''
def prod_maxFR_edges(num_parties, num_edges, num_outcomes):
	return (((num_outcomes**num_edges)**num_edges)**(num_parties-1))**num_parties**num_parties




'''
	For the numerical analysis of the OWB product, this function produces a grid sweep test that indexes
	all the possible configurations of the OWB product within a set of limits, all the while comparing the
	results to the theoretical minimum number of hyperedges achieveable for any compositional product
	that satisfies the No-Signalling condition.
'''
def print_owb_product_grid_sweep(limit_of_parties, limit_of_edges, limit_of_outcomes):
	# For all parametrics
	for parties in range(2,limit_of_parties+1):
		for edges in range(2,limit_of_edges+1):
			for outcomes in range(2,limit_of_outcomes+1):
				print("\n")
				print("Attempting a new configuration:")
				print("Parties: %s | Edges: %s | Outcomes: %s" %( parties, edges, outcomes))
				num_parties = parties
				num_edges = edges
				num_outcomes = outcomes
				y = []
				# Determine the combinations of the parties
				for i in range(num_parties-1, 0, -1):
					y.extend([x for x in itertools.combinations(list(range(num_parties)),i)])
				# Then procure the relations
				total_edges = 0
				for p in y:
					print("\tIndexing relation %s -> %s" % (
						list(p),
						list(filter((None).__ne__,
									[x if x not in list(p) else None for x in list(range(num_parties))]))))
					num_edges_signaller = (num_edges**(len(p)))*((num_outcomes-1)**(len(p)))
					num_edges_signalled = (num_edges**(num_parties-len(p)))-1
					total_add = num_edges_signaller*num_edges_signalled
					print("\t\tAdding %s equalities" % (total_add))
					total_edges += total_add
				# Add on the normalisation edges
				print("\tAdding %s normalisation edges" % ((num_edges**num_parties)))
				print("Total Edges:", (total_edges+(num_edges**num_parties)))
				# Compare to Dr Elie Wolfe's theoretical minimum number of hyperedges (as produced by rref)
				V = num_outcomes
				M = num_edges
				N = num_parties
				print("Elie's Theoretical Minimum:", (((V*M)**N) - (((M*(V-1)) + 1)**N) + 1))
				print("\n")

'''
	This function prints a single test (from within the grid sweep of the OWB product), as a Latex
	expression.
'''
def latex_owb_product(parties, edges, outcomes):
	def ltx_cartesian(arg_input):
		return ' \\times '.join(["E(H_%s)" % (x+1) for x in list(arg_input)])
	def ltx_relation(arg_input):
		return ' \\times '.join(["H_%s" % (x+1) for x in list(arg_input)])
	full_string = ""
	full_string += "\\\\"
	full_string += "\mathrm{Attempting\;a\;new\;configuration\;of\;the\;OWB\;Product}:\\\\"
	full_string += "\;\;\mathrm{For\;a\;configuration\;of\;%s\;parties,\;%s\;hyperedges\;per\;party,\;and\;%s\;outcomes\;per\;hyperedge,}\\\\" %( parties, edges, outcomes)
	num_parties = parties
	num_edges = edges
	num_outcomes = outcomes
	y = []
	for i in range(num_parties-1, 0, -1):
		y.extend([x for x in itertools.combinations(list(range(num_parties)),i)])
	total_edges = 0
	for p in y:
		num_edges_signaller = (num_edges**(len(p)))*((num_outcomes-1)**(len(p)))
		num_edges_signalled = (num_edges**(num_parties-len(p)))-1
		total_add = num_edges_signaller*num_edges_signalled
		total_edges += total_add
		full_string += "\;\;\;\;\mathrm{the\;measurement\;relation\;}\;E_{%s}\\!\\rightarrow\\!E_{%s}\;\mathrm{adds\;%s\;hyperedges,}\\\\" % (ltx_relation(p),ltx_relation(list(filter((None).__ne__,[x if x not in list(p) else None for x in list(range(num_parties))]))),total_add)
	full_string += "\;\;\;\;\mathrm{and\;the\;normalisation\;hyperedges\;that\;form\;the\;Cartesian\;product\;are\;}%s\\\\" % ltx_cartesian([x for x in list(range(num_parties))])
	full_string += "\;\;\;\;\mathrm{which\;require\;%s\;hyperedges}\\\\" % ((num_edges**num_parties))
	full_string += "\mathrm{Total\;hyperedges\;for\;this\;configuration\;is\;%s;}\\\\" % (total_edges+(num_edges**num_parties))
	V = num_outcomes
	M = num_edges
	N = num_parties
	full_string += "\mathrm{Theoretical\;minimum\;is\;}%s" % (((V*M)**N) - (((M*(V-1)) + 1)**N) + 1)
	full_string += "\\\\"
	return full_string

'''
	This function plots the growth of hyperedges within the variants of the FR product, against all known parametrics.
'''
def plot_comparison_of_fr_products():
	fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
	plt.rcParams["font.family"] = 'Times New Roman'
	plt.rcParams["font.weight"] = 'bold'
	plt.rcParams["font.size"] = 24
	plt.rcParams["figure.figsize"] = [24,8]
	################################################
	red_deep = '#872F50'
	red_deft = '#DD4980'
	red_ligt = '#FF679F'
	blu_deft = '#4793CD'
	limparties = 11
	################################################
	minpr = [prod_minFR_edges(x,2,2) for x in range(2,limparties)]
	compr = [prod_commonFR_edges(x,2,2) for x in range(2,limparties)]
	owbpr = [prod_OWB_edges(x,2,2) for x in range(2,limparties)] # The theoretical minimum is also the OWB product
	################################################
	text = "Parties"
	################################################
	ax1.margins(tight=True)
	ax1.grid(which='minor', color='#EEE', linestyle='-', linewidth=0.5)
	ax1.grid(which='major', color='#DDD', linestyle='-', linewidth=1)
	ax1.plot(owbpr, color=red_ligt, marker='D', markersize=7, lw=2, label="Theoretical Minimum")
	ax1.plot(minpr, color=red_deft, marker='D', markersize=7, lw=2, label="Minimal FR Product")
	ax1.plot(compr, color=red_deep, marker='D', markersize=7, lw=2, label="Common Product")
	ax1.set_xticks(range(9))
	ax1.set_xticklabels(range(2,11))
	ax1.set_ylabel("Edges")
	ax1.set_xlabel("%s" % (text))
	ax1.set_yscale('log')
	################################################
	minpr = [prod_minFR_edges(2,x,2) for x in range(2,limparties)]
	compr = [prod_commonFR_edges(2,x,2) for x in range(2,limparties)]
	owbpr = [prod_OWB_edges(2,x,2) for x in range(2,limparties)]
	################################################
	text = "Measurements"
	################################################
	ax2.set_title("Comparison of FR Product Variants Against Theoretical Minimum Number of Hyperedges", y=1.08)
	ax2.margins(tight=True)
	ax2.grid(which='minor', color='#EEE', linestyle='-', linewidth=0.5)
	ax2.grid(which='major', color='#DDD', linestyle='-', linewidth=1)
	ax2.plot(owbpr, color=red_ligt, marker='D', markersize=7, lw=2, label="Theoretical Minimum")
	ax2.plot(minpr, color=red_deft, marker='D', markersize=7, lw=2, label="Minimal FR Product")
	ax2.plot(compr, color=red_deep, marker='D', markersize=7, lw=2, label="Common Product")
	ax2.set_xticks(range(9))
	ax2.set_xticklabels(range(2,11))
	ax2.set_xlabel("%s" % (text))
	ax2.set_yscale('log')
	lgd = ax2.legend(loc=9, bbox_to_anchor=(0.5,-0.3))
	################################################
	minpr = [prod_minFR_edges(2,2,x) for x in range(2,limparties)]
	compr = [prod_commonFR_edges(2,2,x) for x in range(2,limparties)]
	owbpr = [prod_OWB_edges(2,2,x) for x in range(2,limparties)]
	################################################
	text = "Outcomes"
	################################################
	ax3.margins(tight=True)
	ax3.grid(which='minor', color='#EEE', linestyle='-', linewidth=0.5)
	ax3.grid(which='major', color='#DDD', linestyle='-', linewidth=1)
	ax3.plot(owbpr, color=red_ligt, marker='D', markersize=7, lw=2, label="Theoretical Minimum")
	ax3.plot(minpr, color=red_deft, marker='D', markersize=7, lw=2, label="Minimal FR Product")
	ax3.plot(compr, color=red_deep, marker='D', markersize=7, lw=2, label="Common Product")
	ax3.set_xticks(range(9))
	ax3.set_xticklabels(range(2,11))
	ax3.set_xlabel("%s" % (text))
	ax3.set_yscale('log')
	################################################
	fig.subplots_adjust(bottom=0.5)
	plt.draw()
	plt.show()

'''
	This function plots the growth of hyperedges for the minimal FR product and the OWB product, against all known parametrics.
'''
def plot_comparison_of_fr_owb_products():
	fig, (ax1, ax2, ax3) = plt.subplots(1, 3)# plt.figure(figsize=(8,4), dpi=100)
	plt.rcParams["font.family"] = 'Times New Roman'
	plt.rcParams["font.weight"] = 'bold'
	plt.rcParams["font.size"] = 24
	plt.rcParams["figure.figsize"] = [24,8]
	################################################
	red_deep = '#872F50'
	red_deft = '#DD4980'
	red_ligt = '#FF679F'
	blu_deft = '#4793CD'
	limparties = 11
	################################################
	minpr = [prod_minFR_edges(x,2,2) for x in range(2,limparties)]
	owbpr = [prod_OWB_edges(x,2,2) for x in range(2,limparties)]
	text = "Parties"
	################################################
	ax1.margins(tight=True)
	ax1.grid(which='minor', color='#EEE', linestyle='-', linewidth=0.5)
	ax1.grid(which='major', color='#DDD', linestyle='-', linewidth=1)
	ax1.plot(owbpr, color=blu_deft, marker='D', markersize=7, lw=2, label="OWB Product")
	ax1.plot(minpr, color=red_deft, marker='D', markersize=7, lw=2, label="Minimal FR Product")
	ax1.set_xticks(range(9))
	ax1.set_xticklabels(range(2,11))
	ax1.set_ylabel("Edges")
	ax1.set_xlabel("%s" % (text))
	ax1.set_yscale('log')
	################################################
	minpr = [prod_minFR_edges(2,x,2) for x in range(2,limparties)]
	compr = [prod_commonFR_edges(2,x,2) for x in range(2,limparties)]
	owbpr = [prod_OWB_edges(2,x,2) for x in range(2,limparties)]
	################################################
	text = "Measurements"
	################################################
	ax2.set_title("Comparison of Minimal FR Product Against OWB Product Hyperedges", y=1.08)
	ax2.margins(tight=True)
	ax2.grid(which='minor', color='#EEE', linestyle='-', linewidth=0.5)
	ax2.grid(which='major', color='#DDD', linestyle='-', linewidth=1)
	ax2.plot(owbpr, color=blu_deft, marker='D', markersize=7, lw=2, label="OWB Product")
	ax2.plot(minpr, color=red_deft, marker='D', markersize=7, lw=2, label="Minimal FR Product")
	ax2.set_xticks(range(9))
	ax2.set_xticklabels(range(2,11))
	ax2.set_xlabel("%s" % (text))
	ax2.set_yscale('log')
	lgd = ax2.legend(loc=9, bbox_to_anchor=(0.5,-0.3))
	################################################
	minpr = [prod_minFR_edges(2,2,x) for x in range(2,limparties)]
	compr = [prod_commonFR_edges(2,2,x) for x in range(2,limparties)]
	owbpr = [prod_OWB_edges(2,2,x) for x in range(2,limparties)]
	################################################
	text = "Outcomes"
	################################################
	ax3.margins(tight=True)
	ax3.grid(which='minor', color='#EEE', linestyle='-', linewidth=0.5)
	ax3.grid(which='major', color='#DDD', linestyle='-', linewidth=1)
	ax3.plot(owbpr, color=blu_deft, marker='D', markersize=7, lw=2, label="OWB Product")
	ax3.plot(minpr, color=red_deft, marker='D', markersize=7, lw=2, label="Minimal FR Product")
	ax3.set_xticks(range(9))
	ax3.set_xticklabels(range(2,11))
	ax3.set_xlabel("%s" % (text))
	ax3.set_yscale('log')
	################################################
	fig.subplots_adjust(bottom=0.5)
	plt.draw()
	#plt.savefig('all.png', dpi=250,figsize=(8,6))
	plt.show()

'''
	This function produces the number of rows in the row-reduced echelon form matrix of a system consisting
	of three parties, each bearing two measurements with two outcomes. This method of numerical analysis (with
	respect to compositional products of the combinatorial approach) was developed by Sainz and Wolfe in Eq. 2 of "Linear
	Algebra of Hypergraphs", Sainz, Ana Belén, and Elie Wolfe. "Multipartite composition of contextuality
	scenarios." Foundations of Physics 48.8 (2018): 925-953.
'''
def number_of_rows_in_rref_matrix_for_B322():
    # Define the product of two sets of edges in edge-set form
    def product(a,b):
        full_product = []
        for aa in a:
            for bb in b:
                full_product.append(set(['_'.join(sorted(x)) for x in itertools.product(aa,bb)]))
        return full_product
    # Define the relation of two sets of edges in edge-set form
    def relation(a, b):
        def new_edge(w):
            edge = []
            for x in w:
                for y in x[1]:
                    val1 = x[0].split('_')
                    val2 = y.split('_')
                    val1.extend(val2)
                    edge.append('_'.join(sorted(val1)))
            return set(edge)
        full_product = []
        for x in [[a, b]]:
            for a_edge in x[0]:
                for xx in [list(a_edge), list(a_edge)[::-1]]:
                    for i in itertools.product(*([x[1]]*len(xx))):
                        candidate_edge = new_edge(list(zip(xx, i)))
                        if (not candidate_edge in full_product):
                            full_product.append(candidate_edge)
        return full_product
    # Define the edge-sets for three parties of two edges and two measurements each
    a_edges = [['a1','a2'],['a4','a5']]
    b_edges = [['b1','b2'],['b4','b5']]
    c_edges = [['c1','c2'],['c4','c5']]
    # Compute their respective relations
    total_edges = []
    total_edges.extend(relation(product(a_edges,b_edges),c_edges))
    total_edges.extend(relation(product(a_edges,c_edges),b_edges))
    total_edges.extend(relation(product(b_edges,c_edges),a_edges))
    total_edges.extend(relation(c_edges,product(a_edges,b_edges)))
    total_edges.extend(relation(b_edges,product(a_edges,c_edges)))
    total_edges.extend(relation(a_edges,product(b_edges,c_edges)))
    # Procure the columns for their matrix form
    columns = product(product(a_edges,b_edges),c_edges)
    # Procure the rows for the 'No-Signalling' entries of the matrix
    ns_matrix = ''
    for edge_1 in total_edges:
        for edge_2 in total_edges:
            if (edge_1 != edge_2):
                if len(edge_1.intersection(edge_2)) > 0:
                    v_neg = edge_1.difference(edge_2)
                    v_pos = edge_2.difference(edge_1)
                    binary_list = []
                    for c in ([item for sublist in [list(x) for x in columns] for item in sublist]):
                        if (c in list(v_neg)):
                            binary_list.append(1)
                        elif (c in list(v_pos)):
                            binary_list.append(-1)
                        else:
                            binary_list.append(0)
                    ns_matrix += ' '.join([str(x) for x in binary_list])+';'
    ns_matrix = ns_matrix[:-1]
    # Procure the rows for the 'Normalisation' entries of the matrix
    nm_matrix = ''
    for edge_1 in product(product(a_edges,b_edges),c_edges):
        binary_list = []
        for c in ([item for sublist in [list(x) for x in columns] for item in sublist]):
            if (c in edge_1):
                binary_list.append(1)
            else:
                binary_list.append(0)
        nm_matrix += ' '.join([str(x) for x in binary_list])+';'
    nm_matrix = nm_matrix[:-1]
    # Procure the matrix
    all_matrix = np.matrix(ns_matrix+';'+nm_matrix)
    return np.linalg.matrix_rank(all_matrix)