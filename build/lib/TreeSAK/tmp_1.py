import os
import glob
import math
import random
import argparse
import seaborn as sns
from ete3 import Tree
from itolapi import Itol
from PyPDF3.pdf import PageObject
from PyPDF3 import PdfFileWriter, PdfFileReader
from PyPDF3.generic import RectangleObject
# https://pypdf2.readthedocs.io/en/3.0.0/user/adding-pdf-annotations.html

def merge_pdf(pdf_1, pdf_2, margin_size, op_pdf):

    page1 = PdfFileReader(open(pdf_1, "rb"), strict=False).getPage(0)
    page2 = PdfFileReader(open(pdf_2, "rb"), strict=False).getPage(0)

    total_width  = page1.mediaBox.upperRight[0] + page2.mediaBox.upperRight[0] + margin_size*3
    total_height = max([page1.mediaBox.upperRight[1], page2.mediaBox.upperRight[1]]) + margin_size*2
    new_page = PageObject.createBlankPage(None, total_width, total_height)
    new_page.mergeTranslatedPage(page1, margin_size, (total_height-margin_size-page1.mediaBox.upperRight[1]))
    new_page.mergeTranslatedPage(page2, (page1.mediaBox.upperRight[0] + margin_size*2), margin_size)
    output = PdfFileWriter()
    output.addPage(new_page)
    output.write(open(op_pdf, "wb"))


merge_pdf('/Users/songweizhi/Desktop/1.pdf',
          '/Users/songweizhi/Desktop/2.pdf',
          66,
          '/Users/songweizhi/Desktop/merged.pdf')
