#!/usr/bin/env python
# coding: utf-8

# # helpful variables

# In[34]:


standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}

from Bio.Seq import Seq


# # Class Seq

# In[35]:


class seq:
    def __init__(self, name, organism, sequence, type):
        # Call instance attributes
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    # define the info function
    def info(self):
        print(
            self.name + ", " + self.organism + ", " + self.type + ", " + self.sequence
        )
        # print(self.type)
        # print(self.organism)
        # print(self.sequence)

    # define the length function
    def length(self):
        print(len(self.sequence))

    # define the fasta_out function.
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()


# # class protein

# In[92]:


# Write the new protein child class here
class protein(seq):
    def __init__(self, name, organism, sequence, type, size):
        self.size = size
        super().__init__(name, organism, sequence, type)

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "_"
            + self.size
            + "\n"
            + self.sequence
        )
        f.close()

    # To the protein class, add a method called mol_weight, which returns the total molecular
    # weight of the protein sequence. The variable aa_mol_weights in the “helpful variables”
    # file should be helpful. This is a python dictionary of molecular weights for each amino acid
    def mol_weight(self):
        weight_aa = 0
        for i in self.sequence:
            weight_aa = weight_aa + aa_mol_weights[i]
        return weight_aa


# # Class nucleotide

# In[94]:


# Write the new nucleotide class here
class nucleotide(seq):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def gc_content(self):
        countdenom = 0
        countnum = 0
        for i in self.sequence:
            countdenom = countdenom + 1
            if i == "G" or i == "C":
                countnum = countnum + 1
        print((countnum / countdenom) * 100)


# # class DNA


# In[96]:


# Write the DNA class here
class DNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def transcribe(self):
        transcript = self.sequence.replace("T", "U")
        return transcript

    def six_frames(self):
        frame1 = self.sequence
        frame2 = self.sequence[1:]
        frame3 = self.sequence[2:]
        frame4 = self.sequence[::-1]
        frame5 = self.sequence[len(self.sequence) - 1 : 0 : -1]
        frame6 = self.sequence[len(self.sequence) - 2 : 0 : -1]
        print(frame1, frame2, frame3, frame4, frame5, frame6)

    def reverse_complement(self):
        my_DNA = Seq(self.sequence)
        print(my_DNA.reverse_complement())


# ## class RNA

# In[123]:


# Write the RNA class here
class RNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def start(self):
        print(self.sequence.find("AUG"))

    def translate(self):
        startingseq = self.sequence[self.sequence.find("AUG") :]
        returned_amino_acids = []
        for i in range(0, len(startingseq), 3):
            codon = startingseq[i : i + 3]
            print(codon)
            amino_acid = standard_code.get(codon)
            returned_amino_acids.append(amino_acid)
        return returned_amino_acids


# ## Test

# In[125]:


uidA = DNA(
    name="uidA",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
    organism="Bacteria",
    type="DNA",
)


# In[126]:


uidA.fasta_out()


# In[127]:


uidA.six_frames()
uidA.reverse_complement()


# In[128]:


uidA_transcribe = uidA.transcribe()
print(uidA_transcribe)
uid_RNA = RNA(
    name="uidA_RNA", sequence=uidA_transcribe, organism="Bacteria", type="RNA"
)
print(uid_RNA.name, uid_RNA.sequence, uid_RNA.organism, uid_RNA.type)


# In[129]:


uid_RNA.fasta_out()


# In[130]:


uidA_translate = uid_RNA.translate()
uidA_protein = protein(
    name="uidA_protein",
    sequence=uidA_translate,
    organism="Bacteria",
    type="protein",
    size="12",
)
print(uid_RNA.translate())
print(uidA_protein.mol_weight())


# In[131]:


uidA_protein.mol_weight()
