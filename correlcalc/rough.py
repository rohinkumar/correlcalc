            # if value.lower() == 'apdz':
            #         parmetric = APdz
            #         binsparv = binspar
            #         sflag = False
            #         filtermetric = APzdth
            #     elif value.lower() == 'apzdth':
            #         parmetric = APzdth
            #         binsparv = binspar
            #         sflag = False
            #         filtermetric = APzdth
            #     elif value.lower() == 'sflat':
            #         parmetric = flatdistsq
            #         binsparv = binspar**2
            #         filtermetric = flatdistsq
            #         geometry = 'flat'
            #         maxrad = max(binsparv)
            #         # sflag = True
            #     elif value.lower() == 'sopen':
            #         parmetric = opendistsq
            #         binsparv = binspar**2
            #         filtermetric = opendistsq
            #         geometry = 'open'
            #         maxrad = max(binsparv)
            #         # sflag = True
            #     elif value.lower() == 'sclose':
            #         parmetric = closedistsq
            #         binsparv = binspar**2
            #         filtermetric = closedistsq
            #         geometry = 'close'
            #         maxrad = max(binsparv)
            #         # sflag = True
            #     elif value.lower() == 'mu':
            #         parmetric = musq
            #         binsparv = binspar
            #         # sflag = True
            #     elif value.lower() == 'sparf':
            #         parmetric = sparfsq
            #         binsparv = binspar**2
            #         geometry = 'flat'
            #
            #     elif value.lower() == 'sparo':
            #         parmetric = sparosq
            #         binsparv = binspar**2
            #         geometry = 'open'
            #
            #     elif value.lower() == 'sparc':
            #         parmetric = sparcsq
            #         binsparv = binspar**2
            #         geometry = 'close'
            #
            #     else:
            #         print("Incorrect parallel metric argument provided!")
            # elif key.lower() == 'permetric':
            #     if value.lower() == 'apzdth':
            #         permetric = APzdth
            #         binsperv = binsper
            #         maxrad = max(np.sqrt(binsparv**2 + binsperv**2))
            #     elif value.lower() == 'apdz':
            #         permetric = APdz
            #         binsperv = binsper
            #         maxrad = max(np.sqrt(binsparv**2 + binsperv**2))
            #     elif value.lower() == 'sflat':
            #         permetric = flatdistsq
            #         binsperv = binsper**2
            #         filtermetric = flatdistsq
            #         geometry = 'flat'
            #         maxrad = max(binsperv)
            #     elif value.lower() == 'sopen':
            #         permetric = opendistsq
            #         binsperv = binsper**2
            #         filtermetric = opendistsq
            #         geometry = 'open'
            #         maxrad = max(binsperv)
            #     elif value.lower() == 'sclose':
            #         permetric = closedistsq
            #         binsperv = binsper**2
            #         filtermetric = closedistsq
            #         geometry = 'close'
            #         maxrad = max(binsperv)
            #     elif value.lower() == 'mu':
            #         permetric = mu
            #         binsperv = binsper
            #     elif value.lower() == 'sperf':
            #         permetric = sperfsq
            #         binsperv = binsper**2
            #         geometry = 'flat'
            #         maxrad = max(binsparv + binsperv)
            #     elif value.lower() == 'spero':
            #         permetric = sperosq
            #         binsperv = binsper**2
            #         geometry = 'open'
            #         maxrad = max(binsparv + binsperv)
            #     elif value.lower() == 'sperc':
            #         permetric = spercsq
            #         binsperv = binsper**2
            #         geometry = 'close'
            #         maxrad = max(binsparv + binsperv)
            #     else:
            #         print("Incorrect perpendicular metric provided!")
