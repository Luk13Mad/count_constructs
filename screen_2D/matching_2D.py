def find_wMM(template : str, query : str , mismatches_allowed : int) -> int:
    '''searches query in template with specified amount of mismatches allowed.
    Returns -1 if not found, else position of start of first match.
    '''
    if len(query) > len(template):
        return -1
    
    for i,window in sliding_window(template,len(query)):
        if hamming(window,query) <= mismatches_allowed:
            return i
    
    return -1

    
def sliding_window(large_str : str, window_size : int) -> tuple[int,str] :
    '''Sliding window iterator over str.
    Returns window and idx of window start.
    '''
    if len(large_str) == window_size:
        return 0,large_str
    elif len(large_str) < window_size:
        raise ValueError("Window must not be larger than str.")
    
    for i in range(len(large_str) - window_size + 1):
        yield i,large_str[i:i+window_size]


def hamming(str1 : str, str2 : str) -> int:
    '''Hamming distance between str1 and str2'''
    if len(str1) != len(str2):
        raise ValueError(f"Two strings must be same length. str1: {len(str1)} str2: {len(str2)}")
    return sum([x != y for x,y in zip(str1,str2)])
