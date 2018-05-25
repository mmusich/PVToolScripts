import os
import stat
import sys
import time
from os import listdir
from os.path import isfile, join

texTemplate="""
\\documentclass{beamer}
\\usepackage[latin1]{inputenc}
\\usepackage{color}
\\usepackage{longtable}
\\usetheme{default}
\\title[Primary Vertex Validation]{Run by Run Primary Vertex Validation}
\\author{Marco Musich}
\\institute{CP3 Universit\`e catholique de Louvain}
\\date{\\today}

\\begin{document}

\\begin{frame}
\\titlepage
\\end{frame}

\\section{Introduction}
\\begin{frame}{Introduction}
{
\\begin{itemize}
 \\item Run by Run PV Validation
\\end{itemize}
}
\\end{frame}

\\section{Plots}
"""

endTexTemplate="""
\\section{Conclusions}
\\begin{frame}{Conclusions}
{
\\begin{itemize}
\\item Conclusions here
\\end{itemize}
}
\\end{frame} 

\\end{document}
"""

frameTemplate="""
\\begin{frame}{Run number %s}
  \\begin{figure}
    \\centering
    \\includegraphics[width=0.8\\textwidth, height=0.8\\textheight, keepaspectratio=true]{%s}
  \\end{figure}
\\end{frame}
"""

# Script execution.
def main():
    print 'Producing a .tex file from plots...'
    
    mypath="./"

    # Compose .tex frames
    frames = ''
    
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

    ## sort the files
    onlyfiles.sort()
    for ffile in onlyfiles:
        if (ffile.endswith('.png') and not 'L1' in ffile):
            #print ffile
            runnumber = ffile.split('.')[0].split('_')[1]
            frame = frameTemplate % (str(runnumber),str(ffile))
            frames+=frame
            #print frame

    #print frames

    out = texTemplate+frames+endTexTemplate
    print out

    # # Write final .tex file
    ffile = open('presentation.tex', 'w')
    ffile.write(out) 
    # file.write(texTemplate.replace('[frames]', frames).\
    #            replace('[time]', time.ctime()))
    ffile.close()

    # # A script to get from .tex to .pdf
    # pdfScript = open('toPdf.sh', 'w')
    # pdfScript.write(toPdf)
    # pdfScript.close()
    # os.chmod("toPdf.sh", stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)



if __name__ == '__main__':
    main()
