#Basic python test in "hello world" style ...

#function to be tested
def func(x):
    """
    NAME:
        func
    PURPOSE:
        increment input by 1    
    INPUTS:
        x           
    """
    return x+1

# test itself
def test_func():
    #basic assertion
    """
    NAME:
        test_func
    PURPOSE:
        test func    
    INPUTS:
        none           
    """  
    assert func(3) == 4
    


