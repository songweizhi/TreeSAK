from PyPDF3 import PdfFileWriter, PdfFileReader
from PyPDF3.pdf import PageObject


def merge_pdf(pdf_1, pdf_2, op_pdf):

    page1 = PdfFileReader(open(pdf_1, "rb"), strict=False).getPage(0)
    page2 = PdfFileReader(open(pdf_2, "rb"), strict=False).getPage(0)

    total_width = page1.mediaBox.upperRight[0] + page2.mediaBox.upperRight[0]
    total_height = max([page1.mediaBox.upperRight[1], page2.mediaBox.upperRight[1]])

    new_page = PageObject.createBlankPage(None, total_width, total_height)

    # Add first page at the 0,0 position
    new_page.mergePage(page1)
    # Add second page with moving along the axis x
    new_page.mergeTranslatedPage(page2, page1.mediaBox.upperRight[0], 0)

    output = PdfFileWriter()
    output.addPage(new_page)
    output.write(open(op_pdf, "wb"))


pdf_1   = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_hgt_plot_dir/OG0000037_genome_tree_with_HGT.pdf'
pdf_2   = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_hgt_plot_dir/OG0000053_genome_tree_with_HGT.pdf'
pdf_out = '/Users/songweizhi/Desktop/result.pdf'
#merge_pdf(pdf_1, pdf_2, pdf_out)
