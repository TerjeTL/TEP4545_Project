 %%
 %% This is file `asmejour.bst', initially generated with the docstrip utility and extensively customized.
 %%
 %% *** This bibtex style follows ASME journals ***
 %%
 %%
 %% This file is asmejour.bst, fourth 1.04 dated 2020/04/15.
 %%
 %% This file was generated and edited by John Lienhard.
 %% This file is a derivative work under Clause 6 of the LPPL license.
 %%
 %% The generated file and the modifications are based on the formats that ASME is using in
 %% published journal papers in early 2020. Details of the modifications follow.
 %%
 %% Created entry:   @online{...} for web publications (and aliased it to @www and @electronic)
 %% Created fields:  eventdate, urldate, urltype, venue, fourth, fourthtype 
 %% Format for url:  urlprefix = {}, format.web.urldate, bbl.urldate; urltype = {} changes "accessed" to other text. 
 %% Added fields:    venue and eventdate to @inproceedings
 %% Added field:     version={}, format.version to @book, @online, @manual; versiontype = {} changes "Version" to other
 %% Updated doi:     doi url prefix updated to https://doi.org/
 %%
 %% Remove "emphasize" from journal article, and linked doi (or else url) to journal name (or book title) 
 %%   without explicitly showing the urls. The urls may be added to the note = {..} field to display explicitly
 %%   doilink fn is similar to tex.stackexchange.com/questions/76174/set-doi-link-on-journal-name-with-bibtex-bst
 %% Modified fin.entry to omit period after a url, including a separate version for article and book which have
 %%   hyperlinks to titles.
 %%   
 %% Capitalize: Vol., Chap., & No. 
 %% Remove "emphasize" before comma in format.bvolume, 
 %% @incollection:  move chapter and edition, change format.in.ed.booktitle, modify format for "series"
 %%
 %% Modify \eprint to include full url, and provided some macros for archive={..} ... 
 %%   ....this format is uncommon in ASME journals so implemented only as a url (but can easily patch if needed)
 %%
 %% Change format.note capitalizations from "l" to "t" and "u" to "t" respectively.
 %%
 %% Included macros for the titles of ASME's journals, following ASME formats (updated 4/10/20)
 %%
 %%
 %% Modifications in v1.03 (2020/04/12)
 %%
 %% @article:		add warning for empty pages; change ordering of numpages "pages" for eid
 %% @book:          reorder edition, version, and volume 
 %% @inbook:        reorder edition, version, volume, etc.; hyperlink title to doi, url, or eprint 
 %% @incollection:  create format.doi.booktitle to hyperlink title of collection to doi, url, or eprint. 
 %% @inproceedings: move venue, eventdate, and paper number toward end; warn that address should not be used with venue;
 %%					warn for empty venue and empty date;
 %%					create format.pvolume and format.in.ed.proctitle, to permit both [paper] number and volume 
 %%                 print either organization or publisher but not both
 %%					add bbl.paperno = "Paper No." and enable type ={} to change it (this is not backward compatible!)
 %%					hyperlink paper number to doi, or url, or eprint by default.
 %% @inserialcollection: new entry type created to deal with citing articles in edited volumes that  
 %%                      come under different titles in a an ongoing series. This changes the placement
 %%						 of Vol. number in the reference, putting series name and volume in ( ) after title.
 %%						 This calls a new function format.in.ed.serialtitle. Hyperlink title of book.
 %% 					 For most uses, ordinary @incollection will be preferred. 
 %%						 This edit also reverts @incollection and @inproceedings to use Series, Vol. after eds,
 %%						 which is consistent with many ASME references 
 %% @proceedings:   edit order of fields; warn that address should not be used with venue; 
 %%                 print either organization or publisher but not both;
 %%					hyperlink title of proceedings to doi, url, or eprint;
 %%					provide format.alt.link to print the hyperlink if the report number is empty. 
 %% @techrep:       move report number after address; hyperlink report number to doi, url, or eprint;
 %%					provide format.alt.link to print the hyperlink if the report number is empty. 
 %%					provide institution.address function
 %%
 %% Modifications in v1.04 ( )
 %%					add hdl macro as a synonym for hndl
 %%					separate format.serial.volume from format.bvolume; 
 %%					arrange format.bvolume to order vol and series as ", Vol. X of Series Title", which affects
 %%					@book, @inbook, @incollection, and @proceedings.
 %%
 %%
 %%  Author: John H. Lienhard V
 %%          Department of Mechanical Engineering
 %%          Massachusetts Institute of Technology
 %%          Cambridge, MA 02139-4307 USA
 %%
 %========================================================================
 %% LICENSE:
 %%
 %% Copyright (c) 2020 John H. Lienhard
 %%
 %% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
 %% associated documentation files (the "Software"), to deal in the Software without restriction, 
 %% including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 %% and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
 %% subject to the following conditions:
 %%
 %% The above copyright notice and this permission notice shall be included in all copies or 
 %% substantial portions of the Software.
 %%
 %% The software is provided "as is", without warranty of any kind, express or implied, including but 
 %% not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. 
 %% in no event shall the authors or copyright holders be liable for any claim, damages or other liability, 
 %% whether in an action of contract, tort or otherwise, arising from, out of or in connection with the 
 %% software or the use or other dealings in the software.
 %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%
 % The original source files for docstrip were:
 %
 % merlin.mbs  (with options: `seq-no,nm-rev,jnrlst,dt-beg,xmth,tit-qq,atit-u,thtit-a,vol-bf,volp-com,jwdpg,pp-last,num-xser,numser,ser-vol,ser-ed,pre-edn,doi,edby,edbyy,blk-com,in-x,pp,ed,abr,ednx,ord,jabr,revdata,eprint,url,url-blk,eprint,url,url-nt,em-it,nfss,')
 % ----------------------------------------
 % 
 % Copyright 1994-2011 Patrick W Daly
 % ===============================================================
 % IMPORTANT NOTICE:
 % This bibliographic style (bst) file has been generated from one or
 % more master bibliographic style (mbs) files, listed above.
 %
 % This generated file can be redistributed and/or modified under the terms
 % of the LaTeX Project Public License Distributed from CTAN
 % archives in directory macros/latex/base/lppl.txt; either
 % version 1 of the License, or any later version.
 % ===============================================================
 % Name and version information of the main mbs file:
 % \ProvidesFile{merlin.mbs}[2011/11/18 4.33 (PWD, AO, DPC)]
 %   For use with BibTeX version 0.99a or later
 %
 % ===============================================================
 %-------------------------------------------------------------------
 %-------------------------------------------------------------------
 % This bibliography style file is intended for texts in ENGLISH
 % This is a numerical citation style, and as such is standard LaTeX.
 % It requires no extra package to interface to the main text.
 % The form of the \bibitem entries is
 %   \bibitem{key}...
 % Usage of \cite is as follows:
 %   \cite{key} ==>>          [#]
 %   \cite[chap. 2]{key} ==>> [#, chap. 2]
 % where # is a number determined by the ordering in the reference list.
 % The order in the reference list is that by which the works were originally
 %   cited in the text, or that in the database.
 %---------------------------------------------------------------------

ENTRY
  { address
    archive
    author
    booktitle
    chapter
    collaboration
    doi
    edition
    editor
    eid
    eprint
    eventdate
    howpublished
    institution
    journal
    key
    month
    note
    number
    numpages
    organization
    pages
    publisher
    school
    series
    title
    type
    url
    urldate
    urltype
    venue   
    version
    versiontype
    volume
    year
  }
  {}
  { label }
INTEGERS { output.state before.all mid.sentence after.sentence after.block }
FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}
STRINGS { s t}
FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { add.period$ write$
          newline$
          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { add.period$ " " * write$ }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}
FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}
FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {fin.entry}
{ url empty$ 
  { eprint empty$
    {
    add.period$
    write$
    newline$
    }{
    write$
    newline$
    }
    if$
  }{
  write$
  newline$
  }
  if$
}

%% For entries that hyperlink content
FUNCTION {fin.urllink.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}
FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
        'skip$
        { after.sentence 'output.state := }
      if$
    }
  if$
}
FUNCTION {add.blank}
{  " " * before.all 'output.state :=
}

FUNCTION {date.block}
{
  skip$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}
FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}
FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}
FUNCTION {non.stop}
{ duplicate$
   "}" * add.period$
   #-1 #1 substring$ "." =
}

FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}
FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}
FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}
FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}
FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "\textit{" swap$ * "}" * }
  if$
}
FUNCTION {bolden}
{ duplicate$ empty$
    { pop$ "" }
    { "\textbf{" swap$ * "}" * }
  if$
}
FUNCTION {tie.or.space.prefix}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$
}

FUNCTION {capitalize}
{ "u" change.case$ "t" change.case$ }

FUNCTION {space.word}
{ " " swap$ * " " * }
 % Here are the language-specific definitions for explicit words.
 % Each function has a name bbl.xxx where xxx is the English word.
 % The language selected here is ENGLISH
FUNCTION {bbl.and}
{ "and"}

FUNCTION {bbl.etal}
{ "et~al." }

FUNCTION {bbl.editors}
{ "eds." }

FUNCTION {bbl.editor}
{ "ed." }

FUNCTION {bbl.edby}
{ "edited by" }

FUNCTION {bbl.edition}
{ "ed." }

FUNCTION {bbl.volume}
{ "Vol." }

FUNCTION {bbl.of}
{ "of" }

FUNCTION {bbl.number}
{ "No." }

FUNCTION {bbl.nr}
{ "No." }

FUNCTION {bbl.in}
{ "in" }

FUNCTION {bbl.pages}
{ "pp." }

FUNCTION {bbl.page}
{ "p." }

FUNCTION {bbl.eidpp}
{ "pages" }

FUNCTION {bbl.chapter}
{ "Chap." }

FUNCTION {bbl.techrep}
{ "Tech. Rep." }

FUNCTION {bbl.paperno}
{ "Paper No." }

FUNCTION {bbl.mthesis}
{ "Master's thesis" }

FUNCTION {bbl.phdthesis}
{ "Ph.D. thesis" }

FUNCTION {bbl.first}
{ "1st" }

FUNCTION {bbl.second}
{ "2nd" }

FUNCTION {bbl.third}
{ "3rd" }

FUNCTION {bbl.fourth}
{ "4th" }

FUNCTION {bbl.fifth}
{ "5th" }

FUNCTION {bbl.st}
{ "st" }

FUNCTION {bbl.nd}
{ "nd" }

FUNCTION {bbl.rd}
{ "rd" }

FUNCTION {bbl.th}
{ "th" }

FUNCTION {bbl.urldate}
 { "accessed" }
 
FUNCTION {bbl.version}
 { "Version" }

MACRO {jan} {"Jan."}

MACRO {feb} {"Feb."}

MACRO {mar} {"Mar."}

MACRO {apr} {"Apr."}

MACRO {may} {"May"}

MACRO {jun} {"Jun."}

MACRO {jul} {"Jul."}

MACRO {aug} {"Aug."}

MACRO {sep} {"Sep."}

MACRO {oct} {"Oct."}

MACRO {nov} {"Nov."}

MACRO {dec} {"Dec."}

FUNCTION {eng.ord}
{ duplicate$ "1" swap$ *
  #-2 #1 substring$ "1" =
     { bbl.th * }
     { duplicate$ #-1 #1 substring$
       duplicate$ "1" =
         { pop$ bbl.st * }
         { duplicate$ "2" =
             { pop$ bbl.nd * }
             { "3" =
                 { bbl.rd * }
                 { bbl.th * }
               if$
             }
           if$
          }
       if$
     }
   if$
}

MACRO {arxiv} {"https://arxiv.org/abs/"}

MACRO {googlebooks} {"http://books.google.com/books?id="}

MACRO {hdl} {"http://hdl.handle.net/"}
MACRO {hndl} {"http://hdl.handle.net/"}

MACRO {jstor} {"http://www.jstor.org/stable/"}

MACRO {oclc} {"http://www.worldcat.org/oclc/"}

MACRO {pubmed} {"http://www.ncbi.nlm.nih.gov/pubmed/"}


MACRO {amr} {"Appl. Mech. Rev."}

MACRO {lds} {"ASME Letters Dynamic Systems"}

MACRO {jru} {"ASME J. Risk Uncertainty"}

MACRO {jsbc} {"ASME J. Sustainable Bldg. and Cities"}

MACRO {jam} {"ASME J. Appl. Mech."}

MACRO {jbe} {"ASME J. Biomed. Eng."}

MACRO {jcnd} {"ASME J. Comput. Nonlin. Dyn."}

MACRO {jcise} {"ASME J. Comput. Inf. Sci. Eng."}

MACRO {jds} {"J. Dyn. Syst.-T. ASME"}

MACRO {jep} {"ASME J. Electron. Packag."}

MACRO {jert} {"ASME J. Energy Resour."}

MACRO {jegtp} {"ASME J. Eng. Gas Turbines Power"}

MACRO {jem} {"ASME J. Eng. Mater."}

MACRO {jfe} {"ASME J. Fluids Eng."}

MACRO {jfcst} {"ASME J. Fuel Cell Sci. Tech."}

MACRO {jht} {"ASME J. Heat Transfer"}

MACRO {jmse} {"ASME J. Manuf. Sci. E."}

MACRO {jmd} {"ASME J. Mech. Des."}

MACRO {jmr} {"ASME J. Mech. Robot."}

MACRO {jmdv} {"ASME J. Med. Devices"}

MACRO {jmnm} {"ASME J. Micro Nano-Manuf."}

MACRO {jnem} {"ASME J. Nanotech Eng. Med."}

MACRO {joma} {"ASME J. Offshore Mech. Arct."}

MACRO {jpv} {"ASME J. Press. Vess."}

MACRO {jse} {"ASME J. Sol. Energy"}

MACRO {jtsea} {"ASME J. Therm. Sci. Eng. Appl."}

MACRO {jtb} {"ASME J. Tribol."}

MACRO {jtm} {"ASME J. Turbomach."}

MACRO {jva} {"ASME J. Vib. Acoust."}



MACRO {acmcs} {"ACM Comput. Surv."}

MACRO {acta} {"Acta Inf."}

MACRO {cacm} {"Commun. ACM"}

MACRO {ibmjrd} {"IBM J. Res. Dev."}

MACRO {ibmsj} {"IBM Syst.~J."}

MACRO {ieeese} {"IEEE Trans. Software Eng."}

MACRO {ieeetc} {"IEEE Trans. Comput."}

MACRO {ieeetcad}
 {"IEEE Trans. Comput. Aid. Des."}

MACRO {ipl} {"Inf. Process. Lett."}

MACRO {jacm} {"J.~ACM"}

MACRO {jcss} {"J.~Comput. Syst. Sci."}

MACRO {scp} {"Sci. Comput. Program."}

MACRO {sicomp} {"SIAM J. Comput."}

MACRO {tocs} {"ACM Trans. Comput. Syst."}

MACRO {tods} {"ACM Trans. Database Syst."}

MACRO {tog} {"ACM Trans. Graphic."}

MACRO {toms} {"ACM Trans. Math. Software"}

MACRO {toois} {"ACM Trans. Office Inf. Syst."}

MACRO {toplas} {"ACM Trans. Progr. Lang. Syst."}

MACRO {tcs} {"Theor. Comput. Sci."}

FUNCTION {bibinfo.check}
{ swap$
  duplicate$ missing$
    {
      pop$ pop$
      ""
    }
    { duplicate$ empty$
        {
          swap$ pop$
        }
        { swap$
          pop$
        }
      if$
    }
  if$
}
FUNCTION {bibinfo.warn}
{ swap$
  duplicate$ missing$
    {
      swap$ "missing " swap$ * " in " * cite$ * warning$ pop$
      ""
    }
    { duplicate$ empty$
        {
          swap$ "empty " swap$ * " in " * cite$ * warning$
        }
        { swap$
          pop$
        }
      if$
    }
  if$
}
FUNCTION {format.eprint}
{ eprint duplicate$ empty$
    'skip$
    { "\eprint"
      archive empty$
         { "there's an eprint but no archive in " cite$ * warning$ }
        { "[" * archive * "]" * }
      if$
      "{" * swap$ * "}" *
    }
  if$
}

FUNCTION {format.version}
{ version empty$
 {""}
{version "version" bibinfo.check
  versiontype duplicate$ empty$
    { pop$ bbl.version }
    'skip$
  if$
  "versiontype" bibinfo.check
  swap$ duplicate$ empty$
    { pop$ "t" change.case$ }
    { tie.or.space.prefix * * }
  if$
}
if$
}

FUNCTION {format.web.urldate}
{ urldate empty$
  {""}
  {urldate "urldate" bibinfo.check
  urltype duplicate$ empty$
    { pop$ bbl.urldate }
    'skip$
  if$
  "urltype" bibinfo.check
  swap$ duplicate$ empty$
    { pop$ "t" change.case$ }
    { tie.or.space.prefix * * }
  if$
  ", " *
  }
  if$
}

FUNCTION {format.url}
{
  url
  duplicate$ empty$
    { pop$ "" }
    { "\urlprefix\url{" swap$ * "}" * 
    format.web.urldate swap$ *
    }
  if$
}

INTEGERS { nameptr namesleft numnames }


STRINGS  { bibinfo}

FUNCTION {format.names}
{ 'bibinfo :=
  duplicate$ empty$ 'skip$ {
  's :=
  "" 't :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{vv~}{ll}{, f.}{, jj}"
      format.name$
      bibinfo bibinfo.check
      't :=
      nameptr #1 >
        {
          namesleft #1 >
            { ", " * t * }
            {
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              numnames #2 >
                { "," * }
                'skip$
              if$
              t "others" =
                {
                  " " * bbl.etal *
                }
                {
                  bbl.and
                  space.word * t *
                }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
  } if$
}
FUNCTION {format.names.ed}
{
  'bibinfo :=
  duplicate$ empty$ 'skip$ {
  's :=
  "" 't :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{f.~}{vv~}{ll}{, jj}"
      format.name$
      bibinfo bibinfo.check
      't :=
      nameptr #1 >
        {
          namesleft #1 >
            { ", " * t * }
            {
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              numnames #2 >
                { "," * }
                'skip$
              if$
              t "others" =
                {

                  " " * bbl.etal *
                }
                {
                  bbl.and
                  space.word * t *
                }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
  } if$
}
FUNCTION {format.authors}
{ author "author" format.names
    duplicate$ empty$ 'skip$
    { collaboration "collaboration" bibinfo.check
      duplicate$ empty$ 'skip$
        { " (" swap$ * ")" * }
      if$
      *
    }
  if$
}
FUNCTION {get.bbl.editor}
{ editor num.names$ #1 > 'bbl.editors 'bbl.editor if$ }

FUNCTION {format.editors}
{ editor "editor" format.names duplicate$ empty$ 'skip$
    {
      "," *
      " " *
      get.bbl.editor
      *
    }
  if$
}
FUNCTION {format.doi}
{ doi empty$
    { "" }
    {
      "\doi{" doi * "}" *
    }
  if$
}
FUNCTION {doilink}
{ duplicate$ empty$
{ pop$ "" }
{ doi empty$
    { url empty$
        { eprint empty$
        	{ skip$ }
			{  "\href{" archive eprint * * "}{" * swap$ * "}" *  }
	        if$
	    }
        { "\href{" url * "}{" * swap$ * "}" * }
        if$
        }
    { "\href{https://doi.org/" doi * "}{" * swap$ * "}" * }
  if$
}
if$
}

FUNCTION {format.note}
{
 note empty$
    { "" }
    { note #1 #1 substring$
      duplicate$ "{" =
        'skip$
        { output.state mid.sentence =
          { "t" }
          { "t" }
        if$
        change.case$
        }
      if$
      note #2 global.max$ substring$ * "note" bibinfo.check
    }
  if$
}
FUNCTION {format.title}
{ title 
  "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      "\enquote{" swap$ *
      non.stop
        { ",} " * }
        { "} " * }
      if$
    }
  if$
}
FUNCTION {end.quote.title}
{ title empty$
    'skip$
    { before.all 'output.state := }
  if$
}
FUNCTION {output.bibitem}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {n.dashify}
{
  't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
            { "--" *
              t #2 global.max$ substring$ 't :=
            }
            {   { t #1 #1 substring$ "-" = }
                { "-" *
                  t #2 global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}

FUNCTION {word.in}
{ "" }

FUNCTION {format.date}
{
  ""
  duplicate$ empty$
  year  "year"  bibinfo.check duplicate$ empty$
    { swap$ 'skip$
        { "there's a month but no year in " cite$ * warning$ }
      if$
      *
    }
    { swap$ 'skip$
        {
          swap$
          " " * swap$
        }
      if$
      *
    }
  if$
}
FUNCTION {format.btitle}
{ title "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      emphasize
    }
  if$
}
FUNCTION {format.doi.btitle}
{ title "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      doilink emphasize
    }
  if$
}
FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}
FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { bbl.volume volume tie.or.space.prefix "volume" bibinfo.check * *
      series "series" bibinfo.check
      duplicate$ empty$ 'pop$
        { swap$ bbl.of space.word * swap$
          * }
      if$
      "volume and number" number either.or.check
    }
  if$
}
FUNCTION {format.pvolume}
{ volume empty$
    { "" }
    { bbl.volume volume tie.or.space.prefix
      "volume" bibinfo.check * *
      series "series" bibinfo.check
      duplicate$ empty$ 'pop$
        { ", " * swap$ *}
      if$
    }
  if$
}
FUNCTION {format.serial.volume}
{ volume empty$
    { "" }
    { bbl.volume volume tie.or.space.prefix
      "volume" bibinfo.check * *
      series "series" bibinfo.check
      duplicate$ empty$ 'pop$
        { ", " * swap$ *}
      if$
      "volume and number" number either.or.check
    }
  if$
}
FUNCTION {format.number.series}
{ volume empty$
    { number empty$
        { series field.or.null }
        { series empty$
            { number "number" bibinfo.check }
            { output.state mid.sentence =
                { bbl.number }
                { bbl.number capitalize }
              if$
              number tie.or.space.prefix "number" bibinfo.check * *
              bbl.in space.word *
              series "series" bibinfo.check *
            }
          if$
        }
      if$
    }
    { "" }
  if$
}
FUNCTION {is.num}
{ chr.to.int$
  duplicate$ "0" chr.to.int$ < not
  swap$ "9" chr.to.int$ > not and
}

FUNCTION {extract.num}
{ duplicate$ 't :=
  "" 's :=
  { t empty$ not }
  { t #1 #1 substring$
    t #2 global.max$ substring$ 't :=
    duplicate$ is.num
      { s swap$ * 's := }
      { pop$ "" 't := }
    if$
  }
  while$
  s empty$
    'skip$
    { pop$ s }
  if$
}

FUNCTION {convert.edition}
{ extract.num "l" change.case$ 's :=
  s "first" = s "1" = or
    { bbl.first 't := }
    { s "second" = s "2" = or
        { bbl.second 't := }
        { s "third" = s "3" = or
            { bbl.third 't := }
            { s "fourth" = s "4" = or
                { bbl.fourth 't := }
                { s "fifth" = s "5" = or
                    { bbl.fifth 't := }
                    { s #1 #1 substring$ is.num
                        { s eng.ord 't := }
                        { edition 't := }
                      if$
                    }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
  t
}

FUNCTION {format.edition}
{ edition duplicate$ empty$ 'skip$
    {
      convert.edition
      output.state mid.sentence =
        { "l" }
        { "t" }
      if$ change.case$
      "edition" bibinfo.check
      " " * bbl.edition *
    }
  if$
}
INTEGERS { multiresult }
FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
        { #1 'multiresult := }
        { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}
FUNCTION {format.pages}
{ pages duplicate$ empty$ 'skip$
    { duplicate$ multi.page.check
        {
          bbl.pages swap$
          n.dashify
        }
        {
          bbl.page swap$
        }
      if$
      tie.or.space.prefix
      "pages" bibinfo.check
      * *
    }
  if$
}
FUNCTION {format.journal.pages}
{ pages duplicate$ empty$   
       { "empty pages in " cite$ * warning$ }
       'skip$
     if$
     duplicate$ empty$ 'pop$
     { swap$ duplicate$ empty$
        { pop$ pop$ format.pages}
        {
          ", " *
          swap$
          n.dashify
          pages multi.page.check
            'bbl.pages
            'bbl.page
          if$
          swap$ tie.or.space.prefix
          "pages" bibinfo.check
          * *
          *
        }
      if$
    }
  if$
}
FUNCTION {format.journal.eid}
{ eid "eid" bibinfo.check
  duplicate$ empty$ 'pop$
    { swap$ duplicate$ empty$ 'skip$
      {
          ", " *
      }
      if$
      swap$ *
      numpages empty$ 'skip$
        { numpages bbl.eidpp tie.or.space.prefix
          "numpages" bibinfo.check  * *
          " (" swap$ * ")" * *
        }
      if$
    }
  if$
}
FUNCTION {format.vol.num.pages}
{ volume field.or.null
  duplicate$ empty$ 'skip$
    {
      "volume" bibinfo.check
    }
  if$
  bolden
  number "number" bibinfo.check duplicate$ empty$ 'skip$
    {
      swap$ duplicate$ empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
      swap$
      "(" swap$ * ")" *
    }
  if$ *
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    { "" }
    { type empty$
        { bbl.chapter }
        { type "l" change.case$
          "type" bibinfo.check
        }
      if$
      chapter tie.or.space.prefix
      "chapter" bibinfo.check
      * *
    }
  if$
}

FUNCTION {format.booktitle}
{
  booktitle "booktitle" bibinfo.check
  emphasize
}
FUNCTION {format.doi.booktitle}
{
  booktitle "booktitle" bibinfo.check
  doilink emphasize
}
FUNCTION {format.in.ed.booktitle}
{ format.doi.booktitle duplicate$ empty$ 'skip$
    {
      format.edition duplicate$ empty$ 'pop$
        { ", " swap$ * * }
      if$
      editor "editor" format.names.ed duplicate$ empty$ 'pop$
        {
          get.bbl.editor
          swap$ "," *
          " " * swap$ *
          swap$
          "," *
          " " * swap$
          * }
      if$
      format.bvolume duplicate$ empty$ 'pop$
        { ", " swap$ * * }                    
      if$
      word.in swap$ *
    }
  if$
}
FUNCTION {format.in.ed.proctitle}
{ format.booktitle duplicate$ empty$ 'skip$
    {
      format.edition duplicate$ empty$ 'pop$
        { ", " swap$ * * }
      if$
      editor "editor" format.names.ed duplicate$ empty$ 'pop$
        {
          get.bbl.editor
          swap$ "," *
          " " * swap$ *
          swap$
          "," *
          " " * swap$
          * }
      if$
      format.pvolume duplicate$ empty$ 'pop$
        { ", " swap$ * * }                    
      if$
      word.in swap$ *
    }
  if$
}
FUNCTION {format.in.ed.serialtitle}
{ format.doi.booktitle duplicate$ empty$ 'skip$
    {
      format.serial.volume duplicate$ empty$ 'pop$
        { " (" swap$ * * ")" *}                    
      if$
      format.edition duplicate$ empty$ 'pop$
        { ", " swap$ * * }
      if$
      editor "editor" format.names.ed duplicate$ empty$ 'pop$
        {
          get.bbl.editor
          swap$ "," *
          " " * swap$ *
          swap$
          "," *
          " " * swap$
          * }
      if$
      word.in swap$ *
    }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}
FUNCTION {format.thesis.type}
{ type duplicate$ empty$
    'pop$
    { swap$ pop$
      "t" change.case$ "type" bibinfo.check
    }
  if$
}

FUNCTION {format.alt.link}
{ doi empty$
    { url empty$
        { eprint empty$ 
            {   format.note output 
                fin.urllink.entry
            }   
			{   format.note output 
			    format.eprint output
			    fin.entry
			}
	      if$
        }
        { format.note output 
		  format.url output
		  fin.entry
		}
      if$
    }
    { format.doi output 
      format.note output
	  fin.urllink.entry
    }
  if$
}

FUNCTION {format.tr.number}
{ number "number" bibinfo.check
  type duplicate$ empty$
    { pop$ bbl.techrep }
    'skip$
  if$
  "type" bibinfo.check
  swap$ duplicate$ empty$
    { pop$ "t" change.case$ }
    { tie.or.space.prefix * * }
  if$
    duplicate$ empty$ 'skip$
    {
      doilink
    }
  if$
}

FUNCTION {format.paper.number}
{ number "number" bibinfo.check
  type duplicate$ empty$
    { pop$ bbl.paperno }
    'skip$
  if$
  "type" bibinfo.check
  swap$ duplicate$ empty$
    { pop$ "t" change.case$ }
    { tie.or.space.prefix * * }
  if$
    duplicate$ empty$ 'skip$
    {
    doilink
    }
  if$
}
FUNCTION {format.venue}
{ venue "venue" bibinfo.check
  duplicate$ empty$ 'skip$
    { 
    }
  if$
}
FUNCTION {format.eventdate}
{ eventdate "eventdate" bibinfo.check
  duplicate$ empty$ 'skip$
    {
    }
  if$
}

FUNCTION {format.article.crossref}
{
  key duplicate$ empty$
    { pop$
      journal duplicate$ empty$
        { "need key or journal for " cite$ * " to crossref " * crossref * warning$ }
        { "journal" bibinfo.check emphasize word.in swap$ * }
      if$
    }
    { word.in swap$ * " " *}
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  "editor" bibinfo.check
  editor num.names$ duplicate$
  #2 >
    { pop$
      "editor" bibinfo.check
      " " * bbl.etal
      *
    }
    { #2 <
        'skip$
        { editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
            {
              "editor" bibinfo.check
              " " * bbl.etal
              *
            }
            {
             bbl.and space.word
              * editor #2 "{vv~}{ll}" format.name$
              "editor" bibinfo.check
              *
            }
          if$
        }
      if$
    }
  if$
}
FUNCTION {format.book.crossref}
{ volume duplicate$ empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      pop$ word.in
    }
    { bbl.volume
      swap$ tie.or.space.prefix "volume" bibinfo.check * * bbl.of space.word *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { series empty$
            { "need editor, key, or series for " cite$ * " to crossref " *
              crossref * warning$
              "" *
            }
            { series emphasize * }
          if$
        }
        { key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.incoll.inproc.crossref}
{
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { format.booktitle duplicate$ empty$
            { "need editor, key, or booktitle for " cite$ * " to crossref " *
              crossref * warning$
            }
            { word.in swap$ * }
          if$
        }
        { word.in key * " " *}
      if$
    }
    { word.in format.crossref.editor * " " *}
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.org.or.pub}
{ 't :=
  ""
  address empty$ t empty$ and
    'skip$
    {
      t empty$
        { address "address" bibinfo.check *
        }
        { t *
          address empty$
            'skip$
            { ", " * address "address" bibinfo.check * }
          if$
        }
      if$
    }
  if$
}
FUNCTION {format.publisher.address}
{ publisher "publisher" bibinfo.warn format.org.or.pub
}

FUNCTION {format.organization.address}
{ organization "organization" bibinfo.check format.org.or.pub
}

FUNCTION {format.institution.address}
{ institution "institution" bibinfo.check format.org.or.pub
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  format.date "year" output.check
  date.block
  format.title "title" output.check
  end.quote.title
  crossref missing$
    {
      journal
      "journal" bibinfo.check
      doilink
      "journal" output.check
      format.vol.num.pages output
    }
    { format.article.crossref output.nonnull
    }
  if$
  eid empty$
    { format.journal.pages }
    { format.journal.eid }
  if$ 
  format.note output
  fin.urllink.entry
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.date "year" output.check
  date.block
  format.doi.btitle "title" output.check
  crossref missing$
    { format.edition output
      format.version output
      format.bvolume output 
      format.number.series output
      format.publisher.address output
    }
    {
      format.book.crossref output.nonnull
    }
  if$
  format.note output
  fin.urllink.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  format.date output
  date.block
  format.title "title" output.check
  end.quote.title
  howpublished "howpublished" bibinfo.check output
  address "address" bibinfo.check output
  format.doi output
  format.note output
  format.eprint output
  format.url output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.date "year" output.check
  date.block
  format.doi.btitle "title" output.check
  crossref missing$
    {
      format.edition output
      format.version output
      format.bvolume output
      format.number.series output
      format.publisher.address output
      format.chapter.pages "chapter and pages" output.check
    }
    {
      format.chapter.pages "chapter and pages" output.check
      format.book.crossref output.nonnull
    }
  if$
  format.pages "pages" output.check
  format.note output
  fin.urllink.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  format.date "year" output.check
  date.block
  format.title "title" output.check
  end.quote.title
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.number.series output
      format.publisher.address output
      format.chapter.pages output
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  format.pages "pages" output.check
  format.note output
  fin.urllink.entry
}

FUNCTION {inserialcollection}
{ output.bibitem
  format.authors "author" output.check
  format.date "year" output.check
  date.block
  format.title "title" output.check
  end.quote.title
  crossref missing$
    { format.in.ed.serialtitle "booktitle" output.check
      format.number.series output
      format.publisher.address output
      format.chapter.pages output
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  format.pages "pages" output.check
  format.note output
  fin.urllink.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  format.date "year" output.check
  date.block
  format.title "title" output.check
  end.quote.title
  crossref missing$
    { format.in.ed.proctitle "booktitle" output.check
      publisher empty$
        { format.organization.address output }
        { format.publisher.address output }
      if$
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  venue empty$ 
    { venue "venue" output.check}
    { "address and venue" address either.or.check format.venue output }
  if$
  format.eventdate "eventdate" output.check
  number empty$
  	{ format.pages "pages" output.check
	  format.alt.link 
	}
  	{ format.paper.number "number" output.check
  	  format.pages "pages" output.check
      format.note output
      fin.urllink.entry
  	}
  if$
}
FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  author empty$
    { organization "organization" bibinfo.check
      duplicate$ empty$ 'pop$
        { output
          address "address" bibinfo.check output
        }
      if$
    }
    { format.authors output.nonnull }
  if$
  format.date output
  date.block
  format.btitle "title" output.check
  format.edition output
  format.version output
  author empty$
    { organization empty$
        {
          address "address" bibinfo.check output
        }
        'skip$
      if$
    }
    {
      organization "organization" bibinfo.check output
      address "address" bibinfo.check output
    }
  if$
  format.doi output
  format.note output
  format.eprint output
  format.url output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  format.date "year" output.check
  date.block
  format.title
  "title" output.check
  end.quote.title
  bbl.mthesis format.thesis.type output.nonnull
  school "school" bibinfo.warn output
  address "address" bibinfo.check output
  format.doi output
  format.note output
  format.eprint output
  format.url output
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  format.date output
  format.title output
  end.quote.title
  howpublished "howpublished" bibinfo.check output
  format.doi output
  format.note output
  format.eprint output
  format.url output
  fin.entry
  empty.misc.check
}

FUNCTION {online}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.date "year" output.check
  date.block
  format.title
  "title" output.check
  end.quote.title
  format.version output
  organization "organization" bibinfo.warn output
  address "address" bibinfo.check output  
  format.doi output
  format.note output
  format.eprint output
  format.url "url" bibinfo.warn output
  fin.entry
}
FUNCTION {electronic} { online }
FUNCTION {www} { online }

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  format.date "year" output.check
  date.block
  format.title
  "title" output.check
  end.quote.title
  bbl.phdthesis format.thesis.type output.nonnull
  school "school" bibinfo.warn output
  address "address" bibinfo.check output
  format.doi output
  format.note output
  format.eprint output
  format.url output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization "organization" bibinfo.check output
    }
    { format.editors output.nonnull }
  if$
  format.date "year" output.check
  date.block
  format.doi.btitle "title" output.check
  format.bvolume output
  editor empty$
    { publisher empty$
        {  format.number.series output }
        {
          format.number.series output
          format.publisher.address output
        }
      if$
    }
    { publisher empty$
        {
          format.number.series output
          format.organization.address output }
        {
          format.number.series output
          format.publisher.address output
        }
      if$
     }
  if$
  "address and venue" address either.or.check format.venue "venue" output.check
  format.eventdate "eventdate" output.check
  format.note output
  fin.urllink.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  format.date "year" output.check
  date.block
  format.title
  "title" output.check
  end.quote.title
  format.institution.address output
  number empty$
  	{ format.alt.link }
  	{ format.tr.number output.nonnull
  	  format.note output
  	  fin.urllink.entry
  	}
  if$
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  format.date output
  date.block
  format.title "title" output.check
  end.quote.title
  format.doi output
  format.note "note" output.check
  format.eprint output
  format.url output
  fin.entry
}

FUNCTION {default.type} { misc }
READ
STRINGS { longest.label }
INTEGERS { number.label longest.label.width }
FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #1 'number.label :=
  #0 'longest.label.width :=
}
FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
  number.label #1 + 'number.label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
}
EXECUTE {initialize.longest.label}
ITERATE {longest.label.pass}
FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{"  longest.label  * "}" *
  write$ newline$
  "\newcommand{\enquote}[1]{``#1''}"
  write$ newline$
  "\providecommand{\url}[1]{\texttt{#1}}"
  write$ newline$
  "\providecommand{\urlprefix}{}"
  write$ newline$
  "\expandafter\ifx\csname urlstyle\endcsname\relax"
  write$ newline$
  "  \providecommand{\doi}[1]{doi:\discretionary{}{}{}#1}\else"
  write$ newline$
  "  \providecommand{\doi}{doi:\discretionary{}{}{}\begingroup \urlstyle{rm}\Url}\fi"
  write$ newline$
  "\providecommand{\eprint}[2][]{\href{#1#2}{#1#2}}"
  write$ newline$
}
EXECUTE {begin.bib}
EXECUTE {init.state.consts}
ITERATE {call.type$}
FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}
EXECUTE {end.bib}
%% End of customized bst file
%%
%% End of file `asmejour.bst'.
