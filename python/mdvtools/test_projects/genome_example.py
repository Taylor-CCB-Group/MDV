"""
This is being initially written to help me (Peter Todd) better understand the use of genome browser functionality in MDV.
It may serve as an example project, and a reference for regression-testing etc.
That doesn't necessarily mean that the way it is written represents a recommended best-practice - but there is some hope
that at some point it should be.
"""

from mdvtools.conversions import convert_vcf_to_df, convert_vcf_to_mdv
import os

project_folder = os.path.expanduser("~/mdv/genome_browser")

# this is a very trivial example, may want to expand example to download something more interesting
# result will be an MDVProject with a single datasource, with added genome browser capability
vcf_filename = os.path.join(os.path.dirname(__file__), "sample.vcf")

p = convert_vcf_to_mdv(project_folder, vcf_filename)
p.set_editable(True)

df = convert_vcf_to_df(vcf_filename)
p.add_datasource("sample2", df)

link = {"access_data": True}
# p.insert_link ?

p.serve(
    port=5053
)  # as of this writing there was something else hanging out on the default port...
