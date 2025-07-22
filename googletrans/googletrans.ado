program define googletrans

    version 16.1

    syntax [anything], [src(string) dest(string) LISTLanguages]

    if "`listlanguages'" != "" {
        python: list_languages()
        exit
    }
    if `"`src'"' == "" {
        local src "auto"
    }
    if `"`dest'"' == "" {
        local dest "en"
    }
    python: trans(src_text=Macro.getLocal("anything"), src="`src'", dest="`dest'")
end

version 16.1

python:
import asyncio

try:
    from googletrans import Translator, LANGUAGES
except ImportError:
    raise ImportError("you have to install googletrans package first")

from sfi import Data, Macro, ValueLabel, Characteristic, SFIToolkit

def trans(src_text, src="auto", dest="en"):
    translator = Translator()

    if src_text:
        res = asyncio.run(translate_string(src_text, src, dest))
        SFIToolkit.display(res)
        return

    SFIToolkit.displayln("translating variable labels:")
    for vi in range(0, Data.getVarCount()):
        lang = Characteristic.getVariableChar(vi, "lang")
        SFIToolkit.display("  {0} ... ".format(Data.getVarName(vi)))
        if lang != dest:
            res = asyncio.run(translate_string(Data.getVarLabel(vi), src, dest))
            Data.setVarLabel(vi, res.text)
            Characteristic.setVariableChar(vi, "lang", dest)
            SFIToolkit.displayln("done")
        else:
            SFIToolkit.displayln("already in {0}".format(lang))

    SFIToolkit.displayln("translating value labels:")
    for li in ValueLabel.getNames():
        SFIToolkit.display("  {0} ... ".format(li))
        charname = "{0}_lang".format(li)
        lang = Characteristic.getDtaChar(charname)
        if lang != dest:
            try:
                lbls = ValueLabel.getLabels(li)
            except SystemError:
                SFIToolkit.displayln(f"cannot ready label values")
                continue
            vals = ValueLabel.getValues(li)
            results = asyncio.run(translate_string(lbls, src=src, dest=dest))
            for i in range(0, len(lbls)):
                ValueLabel.setLabelValue(li, vals[i], results[i].text)
            Characteristic.setDtaChar(charname, dest)
            SFIToolkit.displayln("done")
        else:
            SFIToolkit.displayln("already in {0}".format(lang))

async def translate_string(src_text, src="auto", dest="en"):
    async with Translator() as translator:
        res = await translator.translate(text=src_text, src=src, dest=dest)

        return res


def list_languages():
    SFIToolkit.displayln("{hline 6}{c TT}{hline 22}")
    SFIToolkit.displayln("{lalign 5:code} {c |} language")
    SFIToolkit.displayln("{hline 6}{c +}{hline 22}")
    for lng in LANGUAGES.items():
        SFIToolkit.displayln("{{lalign 5:{0}}}".format(lng[0]) + " {c |} " + lng[1])
    SFIToolkit.displayln("{hline 6}{c BT}{hline 22}")
end
