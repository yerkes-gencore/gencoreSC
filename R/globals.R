# fixing `plot.cluster.annot.compare: no visible binding for global variable ‘colSum’``; is there a better solution?
# also fixing `EnsDb2GeneName: no visible global function definition for ‘query’` and `EnsDb2GeneName: no visible global function definition for ‘genes’`
utils::globalVariables(c("colSum", "genes", "query"))
