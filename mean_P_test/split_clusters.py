"""
Split the DES SV file into individual files with specific cluster samples.
"""
import fitsio
fname = "sva1_gold_1.0.2_run_redmapper_v6.3.3_lgt20_catalog.fit"
data,header = fitsio.read(fname,header=True)

print header
print data.shape
rich = data['LAMBDA_CHISQ']
z = data['Z_LAMBDA']
print rich.shape
print rich[:10]
print z[:10]
