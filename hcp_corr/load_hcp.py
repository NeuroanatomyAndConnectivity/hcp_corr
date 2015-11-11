
def t_series(subject = "",
             template = None,
             cnt_files=4,
             hemisphere='LH',
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
        'LH' by default, meaning left hemisphere
        
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
    
    for x in xrange(0, cnt_files):
        img = nb.load(files[x])
        
        # find out the indices of brain nodes in hemisphere of interest
        if hemisphere == 'LH':
            # for the left hemisphere : brainModels[0] 
            hem = 0
        
        N_first = img.header.matrix.mims[1].brainModels[hem].indexOffset
        N_cnt = img.header.matrix.mims[1].brainModels[hem].indexCount

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
    # columns are time-series, rows are brain nodes
    return K
