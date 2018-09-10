###get alt fun to do mutation dropplot

GetaltFun <- function(cancer){
  if(cancer == "KIRP"){
    print('Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Splice_Site, Nonstop_Mutation, 
          Translation_Start_Site, Multi_Hit')
    alter_fun = list(
      background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
      },
      Frame_Shift_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
      },
      Frame_Shift_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
      },
      In_Frame_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Del"], col = NA))
      },
      In_Frame_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Ins"], col = NA))
      }, 
      Missense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Missense_Mutation"], col = NA))
      }, 
      Nonsense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
      }, 
      Splice_Site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Splice_Site"], col = NA))
      },
      Nonstop_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
      }, 
      Translation_Start_Site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Translation_Start_Site"], col = NA))
      }, 
      Multi_Hit = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Multi_Hit"], col = NA))
      }
    )
    return(alter_fun)
  }
  if(cancer == "LAML"){
    print('Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Splice_Site, Nonstop_Mutation, 
          Translation_Start_Site, Multi_Hit')
    alter_fun = list(
      background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
      },
      Frame_Shift_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
      },
      Frame_Shift_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
      },
      In_Frame_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Del"], col = NA))
      },
      In_Frame_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Ins"], col = NA))
      }, 
      Missense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Missense_Mutation"], col = NA))
      }, 
      Nonsense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
      }, 
      Splice_Site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Splice_Site"], col = NA))
      },
      Nonstop_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
      }, 
      Translation_Start_Site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Translation_Start_Site"], col = NA))
      }, 
      Multi_Hit = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Multi_Hit"], col = NA))
      }
    )
    return(alter_fun)
  }
  if(cancer == "LUAD"){
    print('Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Splice_Site, Nonstop_Mutation, 
          Translation_Start_Site, Multi_Hit')
    alter_fun = list(
      background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
      },
      Frame_Shift_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
      },
      Frame_Shift_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
      },
      In_Frame_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Del"], col = NA))
      },
      In_Frame_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Ins"], col = NA))
      }, 
      Missense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Missense_Mutation"], col = NA))
      }, 
      Nonsense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
      }, 
      Splice_Site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Splice_Site"], col = NA))
      },
      Nonstop_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
      }, 
      Translation_Start_Site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Translation_Start_Site"], col = NA))
      }, 
      Multi_Hit = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Multi_Hit"], col = NA))
      }
    )
    return(alter_fun)
  }
  if(cancer == "LIHC"){
    print('Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Splice_Site, Nonstop_Mutation, 
          Translation_Start_Site, Multi_Hit')
    alter_fun = list(
      background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
      },
      Frame_Shift_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
      },
      Frame_Shift_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
      },
      In_Frame_Del = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Del"], col = NA))
      },
      In_Frame_Ins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Ins"], col = NA))
      }, 
      Missense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Missense_Mutation"], col = NA))
      }, 
      Nonsense_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
      }, 
      Splice_Site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Splice_Site"], col = NA))
      },
      Nonstop_Mutation = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
      }, 
      Translation_Start_Site = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Translation_Start_Site"], col = NA))
      }, 
      Multi_Hit = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Multi_Hit"], col = NA))
      }
    )
    return(alter_fun)
  } else {
      print('Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Splice_Site, Nonstop_Mutation, 
            Translation_Start_Site, Multi_Hit')
      alter_fun = list(
        background = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
        },
        Frame_Shift_Del = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
        },
        Frame_Shift_Ins = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
        },
        In_Frame_Del = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Del"], col = NA))
        },
        In_Frame_Ins = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Ins"], col = NA))
        }, 
        Missense_Mutation = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Missense_Mutation"], col = NA))
        }, 
        Nonsense_Mutation = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
        }, 
        Splice_Site = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Splice_Site"], col = NA))
        },
        Nonstop_Mutation = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
        }, 
        Translation_Start_Site = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Translation_Start_Site"], col = NA))
        }, 
        Multi_Hit = function(x, y, w, h) {
          grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Multi_Hit"], col = NA))
        }
      )
      return(alter_fun)
    }
}
