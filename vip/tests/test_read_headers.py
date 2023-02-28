from vip.util.read_sequence import read_headers

def test_read_headers():
    headers = ['test_sequence.fasta', 'NZ_CP065712.1']
    filename = 'vip/tests/datatests/test_sequence.fasta'
    res = read_headers(filename)
    assert res == headers
