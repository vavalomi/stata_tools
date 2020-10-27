{smcl}
{* *! version 1.0  27oct2020}{...}
{viewerjumpto "Syntax" "googletrans##syntax"}{...}
{viewerjumpto "Description" "googletrans##description"}{...}
{viewerjumpto "Options" "googletrans##options"}{...}
{viewerjumpto "Installation" "googletrans##installation"}{...}
{viewerjumpto "Examples" "googletrans##examples"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "googletrans python module" "https://py-googletrans.readthedocs.io"}{...}
{vieweralsosee "Google Translate" "https://translate.google.com"}{...}

{title:Title}

{phang}
{bf:googletrans} {hline 2} Translate strings using Google Translate web API


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:googletrans}
[anythong]
[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt src}}the language of the source text{p_end}
{synopt:{opt dest}}the language to translate the source text into{p_end}
{synopt:{opt listl:anguages}}list all available languages{p_end}

{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:googletrans} translates variable and value labes present in the current dataset using Google Translate web API. This is a Stata wrapper of to {browse "https://pypi.org/project/googletrans":googletrans} package in Python. if {cmd:anything} is specified, only that string will be translated. Otherwise, all variable and value labels will be translated and replaced with the translation.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt src} The language of the source text. The value should be one of the language codes from {cmd:listlanguages} option. If a language is not specified, the system will attempt to identify the source language automatically.

{phang}
{opt dest} The language to translate the source text into. The value should be one of the language codes {cmd:listlanguages} option. If not specified, {it:en} will be assumed.

{phang}
{opt listlanguages} Display the list of all available language codes and exit. All other parameters are ignored. 


{marker installation}{...}
{title:Installation}

{pstd}
The command relies on Python integration thus Stata 16.1 or newer is required.

Before you're able to use the comand googletrans python module must be also installed:
{phang}{cmd:> pip install googletrans}{p_end}

If the following command produces no error you have everything set up correctly:
{phang}{cmd:. python which googletrans}{p_end}

{marker examples}{...}
{title:Examples}

{phang}{cmd:. googletrans "გამარჯობა"}{p_end}

{phang}{cmd:. googletrans, dest("es")}{p_end}
