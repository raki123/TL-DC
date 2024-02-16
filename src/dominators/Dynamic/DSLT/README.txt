To compile the code execute "make" via the command line. To run the code for dbs (or dsnca) execute:  ./dbs <input file> <number of runs>

The input file must be in the following format (DIMACS):
• Must contain the header line: 
  p <number of vertices (the vertices must have IDs at most equal to this number)> <number of edges (including edge additions and deletions)> <ID of root vertex> <dummy vertex (ignored in the execution)>
• The edges of the initial graph have the format: 
  a <ID of source vertex> <ID of target vertex>
The edge insertions/deletions (might be intermixed) have the format:
  i <ID of source vertex> <ID of target vertex> (for edge insertions)
  d <ID of source vertex> <ID of target vertex> (for edge deletions)

During the execution the dominator tree is maintained in the array named "idom", where in the position i of the array is the ID of the dominator of the vertex with ID i.
