import numpy as np


def print_table(n,m,row_name,col_name,row_labels,col_labels,A):
   print("\\begin{table}")
   print("\\begin{tabular}[b]{cc|"+("r"*m)+"|}")
   print("\\cline{{3-{:d}}}".format(2+m))
   print("\\multicolumn{1}{l}{} & \\multicolumn{1}{l|}{}"+" & \\multicolumn{{{:d}}}{{c|}}{{{:s}}}\\\\".format(m,col_name)+" \\cline{{3-{:d}}}".format(2+m))
   s = "\\multicolumn{1}{l}{} & \\multicolumn{1}{l|}{} "
   for j in range(m):
      s = s + "& \\multicolumn{{1}}{{c|}}{{{}}} ".format(col_labels[j])
   s = s + "\\\\ \\hline"
   print(s)
   i = 0
   s = "\\multicolumn{{1}}{{|c|}}{{\\multirow{{{:d}}}{{*}}{{\\rotatebox[origin=c]{{90}}{{{:s}}}}}}}".format(n,row_name)
   s = s + " & {}".format(row_labels[i])
   for i in range(m-1):
      s = s + " & \\multicolumn{{1}}{{r|}}{{{}}}".format(A[i,j])
   s = s + " & {} \\\\ \\cline{{2-{}}}".format(A[i,m-1],2+m)
   print(s)
   for i in range(n-1):
      s = "\\multicolumn{{1}}{{|c|}}{{}} & {}".format(row_labels[i+1])
      for j in range(m-1):
         s = s + " &  \\multicolumn{{1}}{{r|}}{{{}}}".format(A[i+1,j])
      s = s + " & {} \\\\".format(A[i+1,m-1])
      if (i==n-2):
         s = s + " \\hline"
      else:
         s = s + " \\cline{{2-{}}}".format(2+m)
      print(s)
      #  \multicolumn{1}{|c|}{} & r2 & \multicolumn{1}{r|}{a21} & \multicolumn{1}{r|}{a22} & \multicolumn{1}{r|}{a23} & a24 \\ \cline{2-6}
   print("\\end{tabular}")
   print("\\end{table}")