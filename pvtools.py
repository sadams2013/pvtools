# Import standard libraries.
import json

# Import external libraries.
import numpy as np
import pandas as pd

class LookupTable:
    """Store liftover data for a gene.

    Parameters
    ----------
    ng : Sequence
        Sequence object for RefSeqGene.
    g7 : Sequence
        Sequence object for GRCh37.
    g8 : Sequence
        Sequence object for GRCh38.
    """
    def __init__(self, ng, g7, g8):
        self.ng = ng
        self.g7 = g7
        self.g8 = g8
        self.df = self._build_lookup_table(ng, g7, g8)

    def _build_lookup_table(self, ng, g7, g8):
        ng_pos1 = np.arange(1, len(ng.seq)+1)
        ng_pos2 = ng_pos1 - ng.data['CDSStarts'][0]
        ng_pos3 = ng.liftover()
        g7_pos = list(range(g7.data['Start'], g7.data['End']+1))
        g8_pos = list(range(g8.data['Start'], g8.data['End']+1))
        allele = np.array(list(ng.seq))
        annot1 = ng.annotate(cds=False)
        annot2 = ng.annotate(cds=True)
        d = {'Start_Position': ng_pos1, 'ATG_Position': ng_pos2,
             'Transcript_Position': ng_pos3, 'GRCh37_Position': g7_pos,
             'GRCh38_Position': g8_pos, 'Allele': allele,
             'Exon_Annotation': annot1, 'CDS_Annotation': annot2}
        return pd.DataFrame(d)

    def to_tsv(self, f):
        self.df.to_csv(f, sep='\t', index=False)

class Sequence:
    """Store data for a single DNA sequence from a FASTA file.

    Parameters
    ----------
    fasta_file : str
        Path to a FASTA file containing the DNA sequence.

    json_file : str
        Path to a JSON file containing metadata for the DNA sequence.
    """
    def __init__(self, fasta_file, json_file=None):
        self.name, self.seq = self._read_fasta_file(fasta_file)
        self.len = len(self.seq)
        self.data = self._read_json_file(json_file)

    def _read_fasta_file(self, fasta_file):
        name = ''
        seq = ''
        with open(fasta_file) as f:
            name = next(f).strip().replace('>', '')
            for line in f:
                seq += line.strip()
        return name, seq

    def _read_json_file(self, json_file):
        if json_file is None:
            return None
        with open(json_file) as f:
            return json.load(f)

    def transcribe(self):
        """Transcribe the DNA sequence.

        Returns
        -------
        str
            mRNA sequence.
        """
        rna = ''
        for i in range(self.data['ExonCount']):
            start = self.data['ExonStarts'][i]
            end = self.data['ExonEnds'][i]
            rna += self.seq[start-1:end]
        return rna

    def get_exon_dataframe(self):
        """Tabulate Exon data.

        Returns
        -------
        pandas.DataFrame
            Dataframe containing Exon data.
        """
        exon_starts = self.data['ExonStarts']
        exon_ends = self.data['ExonEnds']
        exon_names = [f'Exon {x+1}' for x in range(len(exon_starts))]
        intron_starts = [x+1 for x in exon_ends[:-1]]
        intron_ends = [x-1 for x in exon_starts[1:]]
        intron_names = [f'Intron {x+1}' for x in range(len(intron_starts))]
        upstream_start = 1
        upstream_end = exon_starts[0] - 1
        upstream_name = 'Upstream'
        downstream_start = exon_ends[-1] + 1
        downstream_end = len(self.seq)
        downstream_name = 'Downstream'
        starts = exon_starts + intron_starts + [upstream_start, downstream_start]
        ends = exon_ends + intron_ends + [upstream_end, downstream_end]
        names = exon_names + intron_names + [upstream_name, downstream_name]
        df = pd.DataFrame({'Name': names, 'Start': starts, 'End': ends})
        df = df.sort_values('Start')
        df = df.reset_index(drop=True)
        return df

    def get_cds_dataframe(self):
        """Tabulate CDS data.

        Returns
        -------
        pandas.DataFrame
            Dataframe containing CDS data.
        """
        cds_starts = self.data['CDSStarts']
        cds_ends = self.data['CDSEnds']
        cds_names = [f'CDS {x+1}' for x in range(len(cds_starts))]

        intron_starts = [x+1 for x in cds_ends[:-1]]
        intron_ends = [x-1 for x in cds_starts[1:]]
        intron_names = [f'Intron {x+1}' for x in range(len(intron_starts))]

        exon_df = self.get_exon_dataframe()

        upstream_start = 1
        upstream_end = exon_df[exon_df.Name == 'Upstream'].End.values[0]
        upstream_name = 'Upstream'

        utr5_starts = []
        utr5_ends = []
        atg_pos = self.get_atg_pos()
        i = self.get_atg_exon_index()
        for x in range(self.data['ExonCount']):
            start = self.data['ExonStarts'][x]
            end = self.data['ExonEnds'][x]
            if x < i:
                utr5_starts.append(start)
                utr5_ends.append(end)
            elif x == i:
                utr5_starts.append(start)
                utr5_ends.append(atg_pos-1)
            else:
                break
        utr5_names = [f"5' UTR Exon {x+1}" for x in range(len(utr5_starts))]

        utr5_intron_starts = []
        utr5_intron_ends = []
        for utr5_end in utr5_ends[:-1]:
            utr5_intron_starts.append(utr5_end+1)
        for utr5_start in utr5_starts[1:]:
            utr5_intron_ends.append(utr5_start-1)
        utr5_intron_names = [f"5' UTR Intron {x+1}" for x in range(len(utr5_intron_starts))]

        utr3_starts = []
        utr3_ends = []
        stop_pos = self.get_stop_pos()
        i = self.get_stop_exon_index()
        for x in range(self.data['ExonCount']):
            start = self.data['ExonStarts'][x]
            end = self.data['ExonEnds'][x]
            if x < i:
                pass
            elif x == i:
                utr3_starts.append(stop_pos+1)
                utr3_ends.append(end)
            else:
                utr3_starts.append(start)
                utr3_ends.append(end)
        utr3_names = [f"3' UTR Exon {x+1}" for x in range(len(utr3_starts))]

        utr3_intron_starts = []
        utr3_intron_ends = []
        for utr3_end in utr3_ends[:-1]:
            utr3_intron_starts.append(utr3_end+1)
        for utr3_start in utr3_starts[1:]:
            utr3_intron_ends.append(utr3_start-1)
        utr3_intron_names = [f"3' UTR Intron {x+1}" for x in range(len(utr3_intron_starts))]

        downstream_start = exon_df[exon_df.Name == 'Downstream'].Start.values[0]
        downstream_end = len(self.seq)
        downstream_name = 'Downstream'

        starts = cds_starts + intron_starts + utr5_starts + utr5_intron_starts + utr3_starts + utr3_intron_starts + [upstream_start, downstream_start]
        ends = cds_ends + intron_ends + utr5_ends + utr5_intron_ends + utr3_ends + utr3_intron_ends + [upstream_end, downstream_end]
        names = cds_names + intron_names + utr5_names + utr5_intron_names + utr3_names + utr3_intron_names + [upstream_name, downstream_name]
        df = pd.DataFrame({'Name': names, 'Start': starts, 'End': ends})
        df = df.sort_values('Start')
        df = df.reset_index(drop=True)
        return df

    def annotate(self, cds=False):
        if cds:
            df = self.get_cds_dataframe()
        else:
            df = self.get_exon_dataframe()
        annotations = []
        for i, r in df.iterrows():
            n = r.End - r.Start + 1
            annotations += [r.Name] * n
        return annotations

    def liftover(self):
        cds_df = self.get_cds_dataframe()
        cds_pos = []
        cds_sum = 1
        atg_start = self.data['CDSStarts'][0]
        utr5_exon_offset = -1 * self.get_utr5_exon_len()
        utr3_exon_sum = 1
        for i, r in cds_df.iterrows():
            cds_len = r.End - r.Start + 1
            if r.Name.startswith('CDS'):
                cds_pos += list(range(cds_sum, cds_sum + cds_len))
                cds_sum += cds_len
            elif r.Name.startswith('Intron'):
                cds_pos += [f'{cds_sum-1}+{x}' for x in range(1, cds_len+1)]
            elif r.Name == 'Upstream':
                a = self.get_atg_pos() - self.get_utr5_intron_len()
                cds_pos += [x-a for x in range(1, r.End+1)]
            elif r.Name.startswith("5' UTR Exon"):
                a = r.End - r.Start + 1
                cds_pos += [x for x in range(utr5_exon_offset, utr5_exon_offset+a)]
                utr5_exon_offset += a
            elif r.Name.startswith("5' UTR Intron"):
                cds_pos += [f'{utr5_exon_offset-1}+{x}' for x in range(1, cds_len+1)]
            elif r.Name == 'Downstream':
                a = self.get_utr3_exon_len() + 1
                b = r.End - r.Start + 1
                cds_pos += [f'*{x+a}' for x in range(b)]
            elif r.Name.startswith("3' UTR Exon"):
                a = r.End - r.Start + 1
                cds_pos += [f'*{x}' for x in list(range(utr3_exon_sum, utr3_exon_sum+a))]
                utr3_exon_sum += a
            elif r.Name.startswith("3' UTR Intron"):
                cds_pos += [f'*{utr3_exon_sum-1}+{x}' for x in range(1, cds_len+1)]
            else:
                cds_pos += ['.' for x in range(cds_len)]
        if len(cds_pos) != self.len:
            raise ValueError(f"LiftOver length error: expected {self.len} bp, "
                             f"but generated: {len(cds_pos)} bp")
        return [f'c.{x}' for x in cds_pos]

    def get_atg_pos(self):
        return self.data['CDSStarts'][0]

    def get_atg_exon_index(self):
        exon_starts = self.data['ExonStarts']
        exon_ends = self.data['ExonEnds']
        atg_pos = self.get_atg_pos()
        for i in range(self.data['ExonCount']):
            if exon_starts[i] <= atg_pos <= exon_ends[i]:
                return i

    def get_stop_pos(self):
        return self.data['CDSEnds'][-1]

    def get_stop_exon_index(self):
        exon_starts = self.data['ExonStarts']
        exon_ends = self.data['ExonEnds']
        stop_pos = self.get_stop_pos()
        for i in range(self.data['ExonCount']):
            if exon_starts[i] <= stop_pos <= exon_ends[i]:
                return i

    def get_utr5_intron_len(self):
        df = self.get_cds_dataframe()
        df = df[df.Name.str.contains("5' UTR Intron")]
        return sum(df.End - df.Start + 1)

    def get_utr5_exon_len(self):
        df = self.get_cds_dataframe()
        df = df[df.Name.str.contains("5' UTR Exon")]
        return sum(df.End - df.Start + 1)

    def get_utr3_intron_len(self):
        df = self.get_cds_dataframe()
        df = df[df.Name.str.contains("3' UTR Intron")]
        return sum(df.End - df.Start + 1)

    def get_utr3_exon_len(self):
        df = self.get_cds_dataframe()
        df = df[df.Name.str.contains("3' UTR Exon")]
        return sum(df.End - df.Start + 1)
