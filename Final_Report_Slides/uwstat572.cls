\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{uwstat572}[2011/04/16 UW Stat 572 style v0.1] 
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass[11pt]{article}

\RequirePackage{calc}
\RequirePackage[american]{babel}

\setlength{\paperheight}{11in}
\setlength{\paperwidth}{8.5in}
\addtolength{\voffset}{-1.0in}
\addtolength{\hoffset}{-1.0in}
\setlength{\topmargin}{1in}
\setlength{\oddsidemargin}{1in}
\setlength{\evensidemargin}{1in}
\setlength{\textwidth}{\paperwidth - 2in}
\setlength{\textheight}{\paperheight - 2in}
\setlength{\footskip}{36pt}
\setlength{\marginparsep}{0.5cm}
\setlength{\marginparwidth}{1.5cm}
\setlength{\headheight}{0pt}
\setlength{\headsep}{0pt}


\RequirePackage{fancyhdr}
\pagestyle{fancyplain}
\renewcommand{\headrulewidth}{0pt}
\lhead{}
\chead{}
\rhead{}
\lfoot{}
\cfoot{\thepage}
\rfoot{}
