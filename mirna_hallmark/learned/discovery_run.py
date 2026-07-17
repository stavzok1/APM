from mirna_hallmark.learned import discovery as D
if __name__ == "__main__":
    D.run_all(top=150)   # single-pass scans + single-pass deconv on the strongest 150
