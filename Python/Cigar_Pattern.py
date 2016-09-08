import re

# pattern for CIGAR-String                                                       
cigar_pattern = re.compile(r"\d+[MID]{1}")

def parse_cigar(cigar):
  
  pattern = cigar_pattern.findall(cigar)                                       

  for element in pattern:
    #     cig_count        , cig_symbol
    yield int(element[:-1]), element[-1]

# restructure CIGAR-String to avoid repetitive Edit-Operations
def combine_cigar(cigar):

  cig = ""
  pattern = cigar_pattern.findall(cigar)

  tmp = 0 

  for i in range(1,len(pattern)): 

    ccount = int(pattern[i][:-1]) 
    csymbol = pattern[i][-1] 

    prev_ccount = int(pattern[i - 1][:-1]) 
    prev_csymbol = pattern[i - 1][-1] 

    tmp += prev_ccount 

    if csymbol != prev_csymbol: 

      cig += "%d%s" % (tmp, prev_csymbol) 
      if i < len(pattern): 
        tmp = 0 

    if i == len(pattern) - 1: 
      cig += "%d%s" % (tmp + ccount, csymbol) 

  return cig
