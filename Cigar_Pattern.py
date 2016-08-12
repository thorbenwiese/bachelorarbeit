import re

# pattern for CIGAR-String                                                       
cigar_pattern = re.compile(r"\d+[MID]{1}")

def parse_cigar(cigar):
  
  pattern = cigar_pattern.findall(cigar)                                       

  for element in pattern:
    yield int(element[:-1]),element[-1]

"""
def combine_cigar(cigar):                                                
                                                                                  
  cig = ""                                                                     
  tmp = iteration = 0
                                                                                  
  for ccount,csymbol in parse_cigar(cigar):
                                
    if iteration == 0:
      prev_ccount = ccount
      prev_csymbol = csymbol

    else:
      tmp += prev_ccount                                                         
                                                                                  
      if csymbol != prev_csymbol:                                                
                                                                                  
        cig += "%d%s" % (tmp, prev_csymbol)                                      
        tmp = 0                                                                
                                                                                  
      prev_ccount = ccount
      prev_symbol = csymbol

    iteration += 1

  cig += "%d%s" % (tmp + ccount, csymbol)                                  

  return cig 
"""

#def combine(self, cigar, calc_cigar):
def combine_cigar(cigar):

  cig = ""
  pattern = cigar_pattern.findall(cigar)

  tmp = 0 

  for i in range(1,len(pattern)): 

    ccount = int(pattern[i][:-1]) 
    csymbol = pattern[i][-1] 

    prev_ccount = int(pattern[i-1][:-1]) 
    prev_csymbol = pattern[i-1][-1] 

    tmp += prev_ccount 

    if csymbol != prev_csymbol: 

      cig += "%d%s" % (tmp, prev_csymbol) 
      if i < len(pattern): 
        tmp = 0 

    if i == len(pattern)-1: 
      cig += "%d%s" % (tmp + ccount, csymbol) 

  return cig

# calculate costs for CIGAR-String
# coding: M -> '0', I -> '10', D -> '11'
# ==> M 1 Bit, I/D 2 Bit
def calc_bits(cigar):                                                   
     
  bits = 0

  for element in cigar_pattern.findall(cigar):

    ccount = int(element[:-1])                                                 
    csymbol = element[-1]
    
    if csymbol == 'M':
      bits += ccount

    else:
      bits += ccount * 2
                                                                              
  return bits

