
def t_series(subject = "",
             template = None,
             cnt_files=4,
             hemisphere='LH',
             N_first=None,
             N_cnt=None,
             dtype=None,
             normalize=True):
                
    """Load/write Human Connectome Project (HCP) neuroimaging files via NiBabel 
    module. The HCP data is released in GIFTI format (*nii extention) for 
    almost 500 subjects. This script aims to get concetanation of all 
    time-series for each subject. 
    
    subject : string
         subject = data_path + subject_id. 
         e.g. subject = '/a/documents/connectome/_all/100307'
     
    template : string
        Template is the sketch-name of *.nii files (GIFTI format), it is hard-
        coded as template_flag and template_orig...
        
    cnt_files : int
        Number of *.nii files of interest for a subject. The template above
        has 4 forms in total, therefore cnt_files = 4
   
    hemisphere : string

        # LH: CORTEX_LEFT  >> N_first = 0, N_cnt = 29696 
        # RH: CORTEX_RIGHT  >> N_first = 29696, N_cnt = 29716
        # UL: ACCUMBENS_LEFT >> N_first = 59412, N_cnt = 135
        # UR: ACCUMBENS_RIGHT >> N_first = 59547, N_cnt = 140
        # ML: AMYGDALA_LEFT  >> N_first = 59687, N_cnt = 315
        # MR: AMYGDALA_RIGHT  >> N_first = 60002, N_cnt = 332
        # BS: BRAIN_STEM   >> N_first = 60334, N_cnt = 3472
        # CL: CAUDATE_LEFT  >> N_first = 63806, N_cnt = 728
        # CR: CAUDATE_RIGHT >> N_first = 64534, N_cnt = 755
        # EL: CEREBELLUM_LEFT >> N_first = 65289, N_cnt = 8709  
        # ER: CEREBELLUM_RIGHT  >> N_first = 73998, N_cnt = 9144
        # DL: DIENCEPHALON_VENTRAL_LEFT  >> N_first = 83142, N_cnt = 706 
        # DR: DIENCEPHALON_VENTRAL_RIGHT  >> N_first = 83848, N_cnt = 712
        # HL: HIPPOCAMPUS_LEFT  >> N_first = 84560, N_cnt = 764
        # HR: HIPPOCAMPUS_RIGHT  >> N_first = 85324, N_cnt = 795 
        # PL: PALLIDUM_LEFT  >> N_first = 86119, N_cnt = 297
        # PR: PALLIDUM_RIGHT >> N_first = 86416, N_cnt = 260
        # AL: PUTAMEN_LEFT  >> N_first = 86676, N_cnt = 1060
        # AR: PUTAMEN_RIGHT  >> N_first = 87736, N_cnt = 1010
        # TL: THALAMUS_LEFT  >> N_first = 88746, N_cnt = 1288
        # TR: THALAMUS_RIGHT >> N_first = 90034, N_cnt = 1248
        # full : all of them >> N_first = 0, N_cnt = 91282

    K : output, numpy.ndarray
        Concetanation of time-series matrices obtained from each *.nii file. 
    
    References :
        http://www.humanconnectome.org/
        https://github.com/satra/nibabel/tree/enh/cifti2
        right nibabel version to download:        
        $ git clone --branch enh/cifti2 https://github.com/satra/nibabel.git
    
    """
             
    from glob import glob
    import os
    import numpy as np
    import nibabel as nb

    template_flat = 'rfMRI_REST?_??_Atlas_hp2000_clean.dtseries.nii'
    template_orig = 'MNINonLinear/Results/rfMRI_REST?_??/rfMRI_REST?_??_Atlas_hp2000_clean.dtseries.nii'

    # search files in given and default templates
    files = []
    if template != None:
        files = [val for val in sorted(glob(os.path.join(subject, template)))]
    if len(files) == 0:
        files = [val for val in sorted(glob(os.path.join(subject, template_flat)))]
    if len(files) == 0:
        files = [val for val in sorted(glob(os.path.join(subject, template_orig)))]

    if len(files) < cnt_files:
        raise Exception('Not enough files found!')

    files = files[:cnt_files]
    
    # dictionary for brain structures
    label_index = { 'LH':0, 'RH':1, 'UL':2, 'UR':3, 'ML':4, 'MR':5, 'BS':6,
                 'CL':7, 'CR':8, 'EL':9, 'ER':10, 'DL':11, 'DR':12, 'HL':13,
                 'HR': 14, 'PL':15, 'PR':16, 'AL': 17, 'AR':18, 'TL':19, 
                 'TR':20 }
    
    for x in xrange(0, cnt_files):

        img = nb.load(files[x])
        
        # if beginning and end indices given manually        
        if (N_first != None and N_cnt != None):
            
            single_t_series = img.data[:, N_first:N_first+N_cnt].T

        # if a particular brain structure wanted
        elif hemisphere != 'full':
            
            # find out the indices of brain structure of interest
            hem = label_index[hemisphere]
    
            print "BRAIN STRUCTURE: "            
            print img.header.matrix.mims[1].brainModels[hem].brainStructure
            
            N_first = img.header.matrix.mims[1].brainModels[hem].indexOffset
            N_cnt = img.header.matrix.mims[1].brainModels[hem].indexCount
                
            single_t_series = img.data[:, N_first:N_first+N_cnt].T

        # if all brain nodes wanted
        elif hemisphere == 'full':
            
            N_first = 0
            hem = 1
            N_tmp = img.header.matrix.mims[1].brainModels[hem].indexOffset
            N_cnt = img.header.matrix.mims[1].brainModels[hem].indexCount
            N_cnt += N_tmp
            
            single_t_series = img.data[:, N_first:N_first+N_cnt].T

        # length of time series 
        m = single_t_series.shape[1]
        n = single_t_series.shape[0]        
        
        m_last = m
        n_last = n

        if x == 0:
            # In first loop we initialize matrix K to be filled up and returned
            # By default we are using the same dtype like input file (float32)
            init_dtype = single_t_series.dtype if dtype == None else dtype
            K = np.ndarray(shape=[n,m], dtype=init_dtype, order='F')
        else:
            if  m_last != m:
                print "Warning, %s contains time series of different length" % (subject)
            if  n_last != n:
                print "Warning, %s contains different count of brain nodes" % (subject)
            K.resize([n, K.shape[1]+m])

        # concatenation of (normalized) time-series, column-wise
        if normalize:
            mean_series = single_t_series.mean(axis=1)
            std_series = single_t_series.std(axis=1)
            K[:, -m:] = ((single_t_series.T - mean_series) / std_series).T            

        else:
            K[:, -m:] = single_t_series
        del img
        del single_t_series
        
    return K
