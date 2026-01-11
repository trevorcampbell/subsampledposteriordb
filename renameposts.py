import os
import json

for fn in os.listdir('posteriors'):
  ffn = 'posteriors/'+fn
  ffnn = ffn[:-5] + '_subsampled.json'
  os.rename(ffn, ffnn)
  with open(ffnn, 'rb') as f:
    s = json.load(f)
  s['name'] += '_subsampled'
  s['model_name'] += '_subsampled'
  s['dimensions']['SUBIDX'] = 1
  with open(ffnn, 'w') as f:
    json.dump(s, f)

  

  
