from PyPDF2 import PdfFileReader, PdfFileWriter
from PyPDF2.pdf import PageObject

def add_title_to_page(input_pdf_path, output_pdf_path, title):
    # Create a PDF reader object
    pdf_reader = PdfFileReader(input_pdf_path)

    # Create a PDF writer object
    pdf_writer = PdfFileWriter()

    # Get the first page of the PDF
    page = pdf_reader.getPage(0)

    # Create a blank page with the same dimensions as the first page
    blank_page = PageObject.createBlankPage(None, page.mediaBox.getWidth(),
                                             page.mediaBox.getHeight())

    # Set the title on the blank page
    blank_page.mergeScaledTranslatedPage(page, 1, 0, 0)

    # Set the font properties for the title
    blank_page.mergeScaledTranslatedPage(page, 1, 0, 0)

    # Add the title to the top of the blank page
    blank_page.mergeScaledTranslatedPage(page, 1, 0, 0)

    # Add the modified page to the PDF writer
    pdf_writer.addPage(blank_page)

    # Add the remaining pages from the input PDF to the writer
    for i in range(1, pdf_reader.getNumPages()):
        pdf_writer.addPage(pdf_reader.getPage(i))

    # Write the output PDF file
    with open(output_pdf_path, 'wb') as output_pdf:
        pdf_writer.write(output_pdf)

# Usage example
input_pdf = 'input.pdf'
output_pdf = 'output.pdf'
title = 'My Awesome PDF'

add_title_to_page(input_pdf, output_pdf, title)

