
import sys, os,json 
from gensim.models import Word2Vec
from gensim.models import word2vec
import numpy as np 
from Bio import SeqIO 
import gzip
from gensim.models.keyedvectors import KeyedVectors 

l = []
m = []


def fasta_file():
    n = 3
    with open('uniprot_sprot.fasta','rt') as fasta_file:
       for record in SeqIO.parse(fasta_file, "fasta"):
          seq = record.seq
          l.append(seq)
         
          
fasta_file()


def split_ngrams():
    # Fancy way for splitting in kmers
    for seq in l:
      n = 3
      a, b, c = zip(*[iter(seq)]*n), zip(*[iter(seq[1:])]*n), zip(*[iter(seq[2:])]*n)
      str_ngrams = []
      for ngrams in [a,b,c]:
          x = []
          for ngram in ngrams:
            x.append("".join(ngram))
          str_ngrams.append(x)
      m.append(str_ngrams)

split_ngrams()

#corpus
def create_corpus():
    f = open("out.txt",'w')
    for s in m:
      x = s
      for i in x:
         f.write(" ".join(i)+'\n')
    f.close()            

create_corpus()


#model
corpus = word2vec.Text8Corpus("out.txt")     
model = Word2Vec(corpus, size=100, sg=1, window=10, min_count=1, workers=48)
model.save('word2vec.model') # This will be removed later
model = Word2Vec.load("word2vec.model")
model.train(sentences = corpus,total_examples = 10, epochs = 30,report_delay =1,total_words = model.corpus_count)
model.wv.save_word2vec_format("ngram_corpus.txt")

def get_ngram_vectors(file_path):
        ngram_vectors = {}
        vector_length = None
        with open(file_path) as infile:
            for line in infile:
                line_parts = line.rstrip().split()   
                # skip first line with metadata in word2vec text file format
                if len(line_parts) > 2:     
                    ngram, vector_values = line_parts[0], line_parts[1:]
                    x = [float(item) for item in vector_values]
                    
                    ngram_vectors[ngram] = np.array(x, dtype=np.float32)
        return ngram_vectors
        
    
ngram_vectors = get_ngram_vectors("ngram_corpus.txt")

def normalize(x):
    return x / np.sqrt(np.dot(x, x))

n = 3
size = 100

def to_vecs(seq, ngram_vectors):
        ngrams_seq = split_ngrams()


        protvec = np.zeros(size, dtype=np.float32)
        for index in range(len(seq) + 1 - n):
            ngram = seq[index:index + n]
            if ngram in ngram_vectors:
                ngram_vector = ngram_vectors[ngram]
                protvec += ngram_vector
        return normalize (protvec)



def make_protein_vector_for_uniprot(fasta_file, protein_vector_fname, ngram_vectors):
    with open(fasta_file, 'rt') as fasta_file:
        with open(protein_vector_fname, 'w') as output_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                protein_name = record.name.split('|')[-1]
                protein_vector = to_vecs(record.seq, ngram_vectors)

                output_file.write('{}\t{}\n'.format(protein_name, ' '.join(map(str, protein_vector))))
                

make_protein_vector_for_uniprot("uniprot_sprot.fasta", "protein_vector.txt", ngram_vectors)

  


