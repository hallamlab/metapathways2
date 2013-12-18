

import BlastService
from GridParam import GridParam
from BlastBroker import BlastBroker

if __name__=="__main__":
     fastaFile ="/Users/kishori/MetaPathways_2_0/output/test_sample2/orf_prediction/test_sample2.faa" 
     database = "/Users/kishori/MetaPathways_1_0/blastDB/cog-2007-10-30"
     base_output_folder = "/Users/kishori/MetaPathways_2_0/output"
     blastType = "blastp"
     working_directory = "tmp"

     gridParamArray = []
     grid1 = GridParam('sgeadmin', 'sheol')

     grid2 = GridParam('_www', 'shangri-la', working_directory='/common/scratch/')
     grid3 = GridParam('kishori', 'bugaboo.westgrid.ca')
     grid4 = GridParam('kishori', 'jasper.westgrid.ca')

     gridParamArray.append(grid1)
     gridParamArray.append(grid2)
     gridParamArray.append(grid3)
     gridParamArray.append(grid4)
    
     blastBroker = BlastBroker(base_output_folder, 'test_sample2', fastaFile, database, blastType, batchSize=5, working_directory=working_directory, gridParams=gridParamArray)

     print 'calling'
     blastBroker.submit_jobs()
