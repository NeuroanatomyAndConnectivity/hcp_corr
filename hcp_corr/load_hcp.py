
def t_series(data_path,
             subject = "",
             template = None,
             cnt_files=4,
             N_cnt=None,
             N_first=0,
             subject_path=None,
             dtype=None,
             normalize=True):
                
    """Load/write Human Connectome Project (HCP) neuroimaging files via NiBabel 
    module. The HCP data is released in GIFTI format (*nii extention) for 
    almost 500 subjects. This script aims to get concetanation of all 
    time-series for each subject. 
  
    data_path : string
        Local path for the HCP data
        e.g. data_path = '/a/documents/connectome/_all'
    
    subject : string
         Subject ID numbers of individuals scanned for the HCP. 
         e.g. subject='100307'
     
    template : string
        Template is the sketch-name of *.nii files (GIFTI format) of interest
        for a subject. Each subject directory has many *.nii files, but we are 
        pointing the ones of our interest with this template string. 
        (e.g. template = 'rfMRI_REST?_??_Atlas_hp2000_clean.dtseries.nii' ,
        this will be converted into these forms after loop:
       
        rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii
        rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii
        rfMRI_REST2_LR_Atlas_hp2000_clean.dtseries.nii
        rfMRI_REST2_RL_Atlas_hp2000_clean.dtseries.nii )
        
    cnt_files : int
        Number of *.nii files of interest for a subject. The template above
        has 4 forms in total, therefore cnt_files = 4
   
   subject_path : string
       If template requires additional path attachment before it.
       e.g. subject_path = 'MNINonLinear/Results/rfMRI_REST?_??'
    
    (full address of an *nii data = 'data_path/subject/subject_path/template',
    users are recommended to check this manually!!!)

    N_cnt : int
        User-defined number of brain nodes, specially good for test runs.

    N_first : int
        User-defined first brain node we want to get data from.

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
        files = [val for val in sorted(glob(os.path.join(data_path, subject, template)))]
    if len(files) == 0:
        files = [val for val in sorted(glob(os.path.join(data_path, subject, template_flat)))]
    if len(files) == 0:
        files = [val for val in sorted(glob(os.path.join(data_path, subject, template_orig)))]

    if len(files) < cnt_files:
        raise Exception('Not enough files found!')

    files = files[:cnt_files]
    
    for x in xrange(0, cnt_files):
        img = nb.load(files[x])
        
        # brainModels[2] will include both left and right hemispheres

        # count of brain nodes 
        n = img.header.matrix.mims[1].brainModels[2].indexOffset
        n -= N_first
        
        if N_cnt != None:
            n = min(N_cnt, n)
        single_t_series = img.data[:, N_first:N_first+n].T

        # length of time series 
        m = single_t_series.shape[1]

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
            K.resize([n, K.shape[1] + m])

        # concatenation of (normalized) time-series, column-wise
        if normalize:
            
            mean_series = single_t_series.mean(axis=0)
            std_series = single_t_series.std(axis=0)
            K[:, -m:] = (single_t_series - mean_series) / std_series
    
        else:
            K[:, -m:] = single_t_series
        del img
        del single_t_series
    # columns are time-series, rows are brain nodes
    return K
