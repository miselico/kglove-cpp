
TODO
----
Main
* Add command line arguments parsing option

Graph:
* Add option to add inverse relations to graph
* Add option to add the reversed edged to graph (keeping edge label as is)
* Add option to save weighted graph


KGloVe:
* Make different modes accessible
	* Normal
	* With inverse relations
	* With reversed edges
	* 

Random walks (biased RDF2vec)
* Add generation of random walks and saving to file

C++ style issues
* check constness of iterators wherever possible
* make unsigned when it should  be.
* Currently C style file handling is used. Move to C++ style