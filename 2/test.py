
import numpy as np
import pandas as pd
import read_test
info = read_test.init_reading()

from fastdtw import fastdtw
from scipy.spatial.distance import euclidean

from multiprocessing import Pool
import multiprocessing as mp

def applyParallel_1(dfGrouped, func):
    p = Pool(mp.cpu_count())
    ret_list = p.map(func, [group for name, group in dfGrouped])
    return pd.concat(ret_list, axis=1)

def normalize_flux_Stamatis(series):
    maxs = np.max(series["flux"])
    new_flux = series["flux"]/maxs
    series["flux"] = new_flux
    return series


from tqdm import tqdm

for i_c, (end,start) in enumerate(tqdm(read_test.get_chunks(info, 10000))):
    df = (read_test.read_object_by_index_range(info,start,end))

    df["flux"]  = df["flux"]/13675792.0
    means = df.groupby(["object_id"]).agg({'flux':'mean'})
    means_test = means["flux"]
    print (means_test)
    break;
    if(i_c==0):
        means_test.to_csv("test/pos_flux/means_test.csv", mode='a', index=False, header=True)
    else:
        means_test.to_csv("test/pos_flux/means_test.csv", mode='a', index=False, mode=False)