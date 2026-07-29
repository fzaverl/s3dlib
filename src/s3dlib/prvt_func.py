# Copyright (C) Frank Zaverl, Jr.
# See file LICENSE for license information.

import numpy as np
import s3dlib.surface as s3d

''' 
    The following are generic 'private' utilities 
    used for the surface, line, and cloud classes.
'''

def _extract_vals_from_inputCoor(ponts) :
    """
    Separate (N,4) into (N,3) and (N) arrays,

    or just (N,3) and None.
    Used by __init__ and subclassed classmethods pntsurf
    """
    # update for v_1.3.0
    tempCoor = np.array(ponts,dtype=float)
    vVals = None
    if tempCoor.shape[1] == 3 :
        vCoor = tempCoor
    else :
        vCoor = tempCoor[:,:3]
        vVals = tempCoor[:,3]
    return vCoor,vVals


def _getFunctionName(op,xName=None) :
    if xName is not None : return xName
    # assign function name if not a lambda function
    fobar = lambda x : x
    if op.__name__ != fobar.__name__ :
        name = op.__name__
    else :
        name = None
    return name


class KWprocessor():
    """  Access values from defaults or kwargs. """

    def __init__(self,valDict, methodStr=None, checkInputKeys=True, **kargs) :
        """
        Validate input of a method's kwargs.  Internal use ONLY.

        Parameters
        ----------
        valDict : dictionary of argument names, defaults and valid values.
            key -   keyword argument name
            value - single default value, or list of 2D lists to be validated:
                    If a list of list, the structure is:
                    value[0] - ['default value']    <- single value
                    value[1] - [ a, b, c, ..... ]   <- list of a list of valid values

        methodStr : string, optional, default: None
            Identifier of the method called where a validation error occured
            which is shown in the error message.

        checkInputKeys : bool, default: True
            If True, check if all kargs are in the valDict.
            Unknown keys will produce a ValueError.
            If False, ignore kargs not in valDict.

        kargs :
            The dictionary object containing the input keyword:values to be extracted.

        """
        # FutDev: could use a dictionary for value[1] to validate ranges and types.
        self.valDict = valDict
        self.kargs = kargs
        self.method = methodStr
        if checkInputKeys :
            if len(kargs) == 0    : return
            if len(valDict) == 0  : return
            errMsg = "Keyword, {}, is not recognized. \n List of possible keys are: {}"
            if methodStr is not None : errMsg = '['+methodStr+'] - ' + errMsg
            for key in kargs.keys() : 
                if not key in valDict.keys() :
                    raise ValueError( errMsg.format(key, list(valDict.keys())  ) )
        return    

    def _defVal(self,v) :
        value = self.valDict[v]
        if isinstance(value,list) :
            firstVal = value[0]
            return firstVal[0] if isinstance(firstVal,(list,tuple)) else value
        return value

    def _valueKeywordError(self,keyStr,usrVal) :
        dicVal = self.valDict[keyStr]
        errList = dicVal[1]
        for i in range(2,len(dicVal) ) : errList += dicVal[i] 
        errA = "Keyword, {}, value of {} is not recognized.\n List of possible values are: {}"
        errB = "Keyword, {}, value of '{}' is not recognized.\n List of possible values are: {}" 
        errStr =  errB if isinstance(usrVal,str) else errA      
        if self.method is not None : errStr = '['+self.method+'] - ' + errStr
        return errStr.format(keyStr,str(usrVal),errList)
        
    def getVal(self,keyStr,checkValidity=True)  :
        """
        Return the named argument value.

        Parameters
        ----------
        keyStr : string. Keyword

        checkValidity : bool, optional, default: True
            If True, karg value is checked to be valid
            among the list of valid values.

        Raises
        ------
        ValueError
            if invalid keyword assignment.

        Returns
        -------
            value of keyword, or the default value if the keyStr is not found in kargs.

        """
        def listOfList(v) :
            if not isinstance(v,(list,tuple)) : return False
            return isinstance(v[0],(list,tuple)) 
        kargs = self.kargs
        if len(kargs) == 0 : return self._defVal(keyStr)
        if keyStr not in self.valDict : # only for development purposes.
            print("################ CODE ERROR KW001 #####################")
            return None
        if keyStr in kargs : 
            dicVal = self.valDict[keyStr]
            usrVal = kargs[keyStr]
            if listOfList(dicVal) :
                for vLst in dicVal :
                    if usrVal in set(vLst) : return vLst[0] #<< value is the first.
                if not checkValidity : return usrVal
                else :
                    raise ValueError( self._valueKeywordError(keyStr,usrVal))
            return usrVal
        return self._defVal(keyStr)

    def filter(self):
        """
        Removes defaults from kargs.
        Use implies that more than default kargs are to be used,
        thus 'checkInputKeys' should be set to False at __init__.
        """
        for key in self.valDict:
            try :   del self.kargs[ key ]
            except: pass
        return self.kargs


