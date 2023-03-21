
from mlmodel.compute_ml_features import * 






if __name__ == '__main__':
    virus_directory_path = './test_set/virus_sequences/'
    host_directory_path = './test_set/host_sequences/'

    blastn_path = './test_set/StaphStudy_virusvhosts.tsv'
    spacer_path = './test_set/StaphStudy_virusvspacers_blastn.tsv'

    model_path = './vip/gbrt.pkl'


    test = ComputeFeatures(virus_directory_path, host_directory_path)
    test.add_blastn_files(blastn_path, spacer_path)
    test.do_setup()
    test.run_parallel(6)
    test.save_features('saved_features_test.tsv')

