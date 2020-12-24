class Sequence:
    """Store data for a single DNA sequence from a FASTA file.

    Parameters
    ----------
    fasta_json : dict
        JSON object.
    """
    def __init__(self, fasta_json):
        self.data = fasta_json
        self.name, self.seq = self._parse_fasta(self.data['File'])

    def _parse_fasta(self, fasta_file):
        name = ''
        seq = ''
        with open(fasta_file) as f:
            name = next(f).strip().replace('>', '')
            for line in f:
                seq += line.strip()
        return name, seq

    def transcribe(self):
        """Transcribe the DNA sequence.

        Returns
        -------
        str
            mRNA sequence.
        """
        rna = ''
        for i in range(self.data['mRNACount']):
            start = self.data['mRNAStarts'][i]
            end = self.data['mRNAEnds'][i]
            rna += self.seq[start-1:end]
        return rna

    def liftover(self, ng_pos, cds_start):
        result = 0
        if ng_pos < self.data['mRNAStarts'][0]:
            result = ng_pos - self.data['mRNAStarts'][0] - cds_start + 1
        else:
            bp_sum = 0
            for i in range(self.data['CDSCount']):
                start = self.data['CDSStarts'][i]
                end = self.data['CDSEnds'][i]
                if start <= ng_pos <= end:
                    result = ng_pos - start + 1 + bp_sum
                else:
                    bp_sum += end - start + 1
        return result
