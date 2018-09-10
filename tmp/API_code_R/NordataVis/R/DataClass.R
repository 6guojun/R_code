 #=======================================================================
 #OGOCArray: the data or Json contain one gene and one cancer information
 #-----------------------------------------------------------------------
 OGOCArray <- setClass("OCOGArray", slots = list(MData = "list", GIList = "list", className = "charecter"),
                       contains = "eSet"
                       )
 