def patterncount(text, pattern):
    count = 0
    i = 0
    end = len(text) - len(pattern)
    for i in range(end+1):
        if text[i:i+len(pattern)] == pattern:
            count +=1
        i +=1
    return count
    
