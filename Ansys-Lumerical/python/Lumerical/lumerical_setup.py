# Copyright (c) 2021, Ansys, Inc. All rights reserved.


import pya


def register_lumerical_menu():
    import os

    menu = pya.Application.instance().main_window().menu()

    s_1 = "lumerical_menu"
    if not(menu.is_menu(s_1)):
        menu.insert_menu("help_menu", s_1, "Ansys Lumerical")

    s_2 = "interconnect"
    if not(menu.is_menu(s_1 + "." + s_2)):
        menu.insert_menu(s_1 + ".end", s_2, "Circuits Simulation")
    
