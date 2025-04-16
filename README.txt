31 July 2020, Version 5.0
==================================
Files in the 
American Meteorological Society
LaTeX Package
==================================
Copyright 2020, American Meteorological Society

This work may be distributed and/or modified under the
conditions of the LaTeX Project Public License, either version 1.3
of this license or (at your option) any later version.
The latest version of this license is in
http://www.latex-project.org/lppl.txt
and version 1.3 or later is part of all distributions of LaTeX
version 2005/12/01 or later.

This work has the LPPL maintenance status `maintained'.
The Current Maintainer of this work is the American Meteorological Society.
This work consists of the files ametsocV5.cls, ametsoc2014.bst, amsdocsV5.pdf, and AMS_RefsV5.pdf.
==================================
You will be provided with a tarred, zipped LaTeX package containing 
17 files.

==================================
Basic Style File: ametsocV5.cls. 

The file ametsocV5.cls is the manuscript style file.  

Use \documentclass{ametsocV5} for your .tex document to
generate a PDF that follows all AMS guidelines for submission and peer review.  

Use \documentclass[twocol]{ametsocV5} for your .tex document
to generate a PDF that closely
follows the layout of an AMS journal page, including single spacing and two
columns.  This journal style PDF is only for the author's personal use, and
any papers submitted in this style will not be accepted.  

Always use \documentclass{ametsocV5} 
when generating a PDF for submission to the AMS. 

==================================
Template:

templateV5.tex is a file for the author to use to create their paper.

The file provides a basic blank template with some
section headings for authors to more easily enter their manuscript.

==================================
Bibliography Files:

ametsoc2014.bst, references.bib, and AMS_RefsV5.pdf  

ametsoc2014.bst is the bibliography style file. 

references.bib should be altered with your own bibliography information. 

AMS_RefsV5.pdf explains AMS reference style and contains detailed examples. 
Reference list included within is obtained from database2020.bib.

==================================
Documentation: 

Found in amsdocsV5.pdf, amspaperV5.tex/amspaperV5.pdf, and in this file, readme.txt.

==================================
Samples:
amssamp1V5.tex
amssamp1V5.pdf
amssamp2V5.tex
amssamp2V5.pdf
database2020.bib
FigOne.pdf
FigTwo.pdf
figure01.pdf
==================================
Help for Authors
==================================

For questions or problems related to submitting your LaTeX manuscript 
to the AMS, see the AMS LaTeX Submission Info webpage:

www.ametsoc.org/pubslatex

or contact latex@ametsoc.org

============================
Version History: 

 Version history

 7 May 2014 -- Nicole Rietmann, AMS
 Updated to version 4.1
 Changed draft line spacing from 1.66 to 2
 Added \bibpunct command to template

 12 May 2014 -- Nicole Rietmann, AMS
 Updated to version 4.2
 Added version and date to footer
 Changed line spacing from 2 to 2.25
 Added if/else statements for fig/table line spacing (\baselinestretch)
 Decreased vertical space after section heads
 Removed en dashes around page numbers
 Removed line after abstract
 Decreased vertical space after appendix title

 16 May 2014 -- Nicole Rietmann, AMS
 Updated to version 4.3
 Added command to fix equation line spacing

 19 May 2014 -- Nicole Rietmann, AMS 
 Updated to version 4.3.1
 Removed unnecessary \usepackage{tikz} command

 25 August 2014 -- Nicole Rietmann, AMS
 Updated to version 4.3.2
 Correction of \BAMS to \bams for journal command
 Addition of JAS to journal list in .cls and template
 Clarified appendix fig/table placement in template and amsdocs.pdf
 Addition of BAMS capsule instructions in template and amsdocs.pdf

 2 January to 31 July 2020 -- Nicole Rietmann, AMS
 Updated to version 5.0
 Addition of commands \usepackage{newtxtext} and \usepackage{newtxmath} to fix bug with some Greek letters not showing up properly in Times font
 Moved math commands to templateV5.tex file for greater visibility and easier manipulation, and added optional hmmax and bmmax commands to prevent "too many math alphabets" error.
 Addition of \usepackage{epstopdf} to resolve issues with certain programs not automatically converting eps files to pdf
 Made abstract margins wider
 Removed separate email command, email now included in corresponding author info
 Removed option to select journal name for two-column
 Added statement for preprint use of two-column version
 Updated references.bib, database2020.bib, and documentation to fix errors and include updates to AMS reference style and formatting
 Added instructions and commands to template and documentation for significance statement and data availability statement
 Removed paragraph indentation at beginning of abstract
