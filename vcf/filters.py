def sample_names_in(sample_names):
    return lambda (sample_name, genotype): sample_name in sample_names

def rsid_in(rsids):
    return lambda (rsid): rsid in rsids
