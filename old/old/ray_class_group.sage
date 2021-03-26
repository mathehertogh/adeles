class RayClassGroup(AbelianGroup_class):
    def __init__(self, number_field, mod_ideal = 1, mod_archimedean = None):
        if mod_archimedean == None:
            mod_archimedean = [0] * len(number_field.real_places())

        bnf = gp(number_field.pari_bnf())
        # Use PARI to compute ray class group
        bnr = bnf.bnrinit([mod_ideal, mod_archimedean],1)
        invariants = bnr[5][2]         # bnr.clgp.cyc
        invariants = [ ZZ(x) for x in invariants ]

        AbelianGroup_class.__init__(self, len(invariants), invariants)
        self.__number_field = number_field
        self.__bnr = bnr

    def __call__(self, *args, **kwargs):
        return group.Group.__call__(self, *args, **kwargs)

    def _element_constructor_(self, *args, **kwargs):
        if isinstance(args[0], AbelianGroupElement):
            return AbelianGroupElement(self, args[0])
        else:
            I = self.__number_field.ideal(*args, **kwargs)

            # Use PARI to compute class of given ideal
            g = self.__bnr.bnrisprincipal(I)[1]
            g = [ ZZ(x) for x in g ]
            return AbelianGroupElement(self, g)