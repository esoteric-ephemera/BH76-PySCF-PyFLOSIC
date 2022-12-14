
def indx_to_name():
    # just a collection of the chemical reactions in BH76 and BH76RC

    bh76_d_odd = {
        1 : 'H + N$_2$O $\\rightarrow$ N$_2$ + OH',
        3 : 'H + HF $\\rightarrow$ HF + H',
        5 : 'H + HCl $\\rightarrow$ HCl + H',
        7 : 'CH$_3$F + H $\\rightarrow$ CH$_3$ + HF',
        9 : 'F$_2$ + H $\\rightarrow$ F + HF',
        11 : 'CH$_3$ + ClF $\\rightarrow$ CH$_3$F + Cl',
        13 : 'CH$_3$F + F$^-$ $\\rightarrow$ FCH$_3$ + F$^-$',
        15 : 'CH$_3$F $\\hdots$ F$^-$ $\\rightarrow$ FCH$_3$ $\\hdots$ F$^-$',
        17 : 'CH$_3$Cl + Cl$^-$ $\\rightarrow$ ClCH$_3$ + Cl$^-$',
        19 : 'CH$_3$Cl $\\hdots$ Cl$^-$ $\\rightarrow$ ClCH$_3$ $\\hdots$ Cl$^-$',
        21 : 'CH$_3$Cl + F$^-$ $\\rightarrow$ CH$_3$F + Cl$^-$',
        23 : 'F$^-$ $\\hdots$ CH$_3$Cl $\\rightarrow$ CH$_3$F $\\hdots$ Cl$^-$',
        25 : 'CH$_3$F + OH$^-$ $\\rightarrow$ HOCH$_3$ + F$^-$',
        27 : 'CH$_3$F $\\hdots$ OH$^-$ $\\rightarrow$ HOCH$_3$ $\\hdots$ F$^-$',
        29 : 'H + N$_2$ $\\rightarrow$ HN$_2$',
        31 : 'H + CO $\\rightarrow$ HCO',
        33 : 'H + C$_2$H$_4$ $\\rightarrow$ CH$_3$CH$_2$',
        35 : 'C$_2$H$_4$ + CH$_3$ $\\rightarrow$ CH$_3$CH$_2$CH$_2$',
        37 : 'HCN $\\rightarrow$ HNC',
        39 : 'H + HCl $\\rightarrow$ H$_2$ + Cl',
        41 : 'OH + H$_2$ $\\rightarrow$ H + H$_2$O',
        43 : 'CH$_3$ + H$_2$ $\\rightarrow$ CH$_4$ + H',
        45 : 'CH$_4$ + OH $\\rightarrow$ CH$_3$ + H$_2$O',
        47 : 'H$_2$ + H $\\rightarrow$ H + H$_2$',
        49 : 'NH$_3$ + OH $\\rightarrow$ NH$_2$ + H$_2$O',
        51 : 'CH$_3$ + HCl $\\rightarrow$ CH$_4$ + Cl',
        53 : 'C$_2$H$_6$ + OH $\\rightarrow$ C$_2$H$_5$ + H$_2$O',
        55 : 'F + H$_2$ $\\rightarrow$ H + HF',
        57 : 'O + CH$_4$ $\\rightarrow$ OH + CH$_3$',
        59 : 'H + PH$_3$ $\\rightarrow$ H$_2$ + PH$_2$',
        61 : 'H + OH $\\rightarrow$ O + H$_2$',
        63 : 'H + H$_2$S $\\rightarrow$ H$_2$ + HS',
        65 : 'O + HCl $\\rightarrow$ OH + Cl',
        67 : 'NH$_2$ + CH$_3$ $\\rightarrow$ NH + CH$_4$',
        69 : 'NH$_2$ + C$_2$H$_5$ $\\rightarrow$ NH + C$_2$H$_6$',
        71 : 'NH$_2$ + C$_2$H$_6$ $\\rightarrow$ NH$_3$ + C$_2$H$_5$',
        73 : 'NH$_2$ + CH$_4$ $\\rightarrow$ NH$_3$ + CH$_3$',
        75 : '\\textit{s-trans cis-}C$_5$H$_8$ $\\rightarrow$ \\textit{s-trans cis-}C$_5$H$_8$'
    }
    bh76_d = {}
    for irx in bh76_d_odd:
        bh76_d[irx] = bh76_d_odd[irx]
        bh76_d[irx+1] = '{:} reverse'.format(irx)

    bh76rc_d = {
        1 : 'H + N$_2$O $\\rightarrow$ OH + N$_2$',
        2 : 'H + CH$_3$F $\\rightarrow$ HF + CH$_3$',
        3 : 'H + F$_2$ $\\rightarrow$ HF + F',
        4 : 'CH$_3$ + ClF $\\rightarrow$ CH$_3$F + Cl',
        5 : 'CH$_3$Cl + F$^-$ $\\rightarrow$ CH$_3$F + Cl$^-$',
        6 : 'F$^-$ $\\hdots$ CH$_3$Cl $\\rightarrow$ CH$_3$F $\\hdots$ Cl$^-$',
        7 : 'CH$_3$F + OH$^-$ $\\rightarrow$ HOCH$_3$ + F$^-$',
        8 : 'CH$_3$F $\\hdots$ OH$^-$ $\\rightarrow$ HOCH$_3$ $\\hdots$ F$^-$',
        9 : 'H + N$_2$ $\\rightarrow$ HN$_2$',
        10 : 'H + CO $\\rightarrow$ HCO',
        11 : 'H + C$_2$H$_4$ $\\rightarrow$ C$_2$H$_5$',
        12 : 'C$_2$H$_4$ + CH$_3$ $\\rightarrow$ CH$_3$CH$_2$CH$_2$',
        13 : 'HCN $\\rightarrow$ HNC',
        14 : 'H + HCl $\\rightarrow$ H$_2$ + Cl',
        15 : 'H$_2$ + OH $\\rightarrow$ H + H$_2$O',
        16 : 'H$_2$ + CH$_3$ $\\rightarrow$ H + CH$_4$',
        17 : 'OH + CH$_4$ $\\rightarrow$ H$_2$O + CH$_3$',
        18 : 'OH + NH$_3$ $\\rightarrow$ H$_2$O + NH$_2$',
        19 : 'CH$_3$ + HCl $\\rightarrow$ CH$_4$ + Cl',
        20 : 'OH + C$_2$H$_6$ $\\rightarrow$ H$_2$O + C$_2$H$_5$',
        21 : 'H$_2$ + F $\\rightarrow$ H + HF',
        22 : 'O + CH$_4$ $\\rightarrow$ OH + CH$_3$',
        23 : 'H + PH$_3$ $\\rightarrow$ H$_2$ + PH$_2$',
        24 : 'H + OH $\\rightarrow$ H$_2$ + O',
        25 : 'H + H$_2$S $\\rightarrow$ H$_2$ + HS',
        26 : 'O + HCl $\\rightarrow$ OH + Cl',
        27 : 'NH$_2$ + CH$_3$ $\\rightarrow$ NH + CH$_4$',
        28 : 'NH$_2$ + C$_2$H$_5$ $\\rightarrow$ NH + C$_2$H$_6$',
        29 : 'NH$_2$ + C$_2$H$_6$ $\\rightarrow$ NH$_3$ + C$_2$H$_5$',
        30 : 'NH$_2$ + CH$_4$ $\\rightarrow$ NH$_3$ + CH$_3$'
    }

    return bh76_d, bh76rc_d
