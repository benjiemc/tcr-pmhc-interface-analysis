'''Constants and functions for annotating sequences as CDR domains in T cell receptors.'''
IMGT_CDR1 = set(range(27, 38 + 1))
'''IMGT residue numbers corresponding to CDR 1 domains.'''
IMGT_CDR2 = set(range(56, 65 + 1))
'''IMGT residue numbers corresponding to CDR 2 domains.'''
IMGT_CDR3 = set(range(105, 117 + 1))
'''IMGT residue numbers corresponding to CDR 3 domains.'''
IMGT_CDR = IMGT_CDR1.union(IMGT_CDR2).union(IMGT_CDR3)
'''IMGT residue numbers corresponding to all CDR domains.'''

IMGT_VARIABLE_DOMAIN = set(range(0, 128 + 1))
'''Variable domain range for IMGT numbered TCR structures.'''


def assign_cdr_number(imgt_id: str | int | None) -> int | None:
    '''
    Map imgt_id to CDR domains, return number associated with domain or return None if input is not in a CDR
    domain.
    '''
    if not imgt_id:
        return None

    if type(imgt_id) is not int:
        seq_id = int(''.join([char for char in imgt_id if char.isnumeric()]))

    else:
        seq_id = imgt_id

    if seq_id in IMGT_CDR1:
        return 1

    if seq_id in IMGT_CDR2:
        return 2

    if seq_id in IMGT_CDR3:
        return 3

    return None
