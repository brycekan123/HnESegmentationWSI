import numpy as np
import cv2
from matplotlib import pyplot as plt

def norm_HnE(img, Io=240, alpha=1, beta=0.15,process=True):
    # Reference H&E OD matrix
    HERef = np.array([[0.5626, 0.2159],
                      [0.7201, 0.8012],
                      [0.4062, 0.5581]])
    
    # Reference maximum stain concentrations for H&E
    maxCRef = np.array([1.9705, 1.0308])
    
    # Extract dimensions
    h, w, c = img.shape
    
    # Reshape image
    img_reshaped = img.reshape((-1, 3))
    
    # Calculate optical density
    OD = -np.log10((img_reshaped.astype(float) + 1) / Io)
    
    # Remove transparent pixels
    ODhat = OD[~np.any(OD < beta, axis=1)]
    
    # Check if we have enough tissue pixels for stain separation
    if len(ODhat) < 100:  # Minimum threshold
        print(f"Warning: Only {len(ODhat)} tissue pixels found. Skipping stain separation.")
        process = False
        return img, img, img,process  # Return original image for all three outputs
    
    try:
        # Calculate SVD on the OD tuples
        eigvals, eigvecs = np.linalg.eigh(np.cov(ODhat.T))
        
        # Project on the plane spanned by the eigenvectors corresponding to the two largest eigenvalues
        That = ODhat.dot(eigvecs[:,1:3])
        
        # Calculate angle of each point wrt the first SVD direction
        phi = np.arctan2(That[:,1],That[:,0])
        
        minPhi = np.percentile(phi, alpha)
        maxPhi = np.percentile(phi, 100-alpha)
        
        vMin = eigvecs[:,1:3].dot(np.array([(np.cos(minPhi), np.sin(minPhi))]).T)
        vMax = eigvecs[:,1:3].dot(np.array([(np.cos(maxPhi), np.sin(maxPhi))]).T)
        
        # Heuristic to make hematoxylin first and eosin second
        if vMin[0] > vMax[0]:    
            HE = np.array((vMin[:,0], vMax[:,0])).T
        else:
            HE = np.array((vMax[:,0], vMin[:,0])).T
        
        # Determine concentrations of the individual stains
        Y = np.reshape(OD, (-1, 3)).T
        C = np.linalg.lstsq(HE, Y, rcond=None)[0]
        
        # Normalize stain concentrations
        maxC = np.array([np.percentile(C[0,:], 99), np.percentile(C[1,:], 99)])
        tmp = np.divide(maxC, maxCRef)
        C2 = np.divide(C, tmp[:, np.newaxis])
        
        # Recreate the normalized image using reference mixing matrix 
        Inorm = np.multiply(Io, np.exp(-HERef.dot(C2)))
        Inorm[Inorm>255] = 254
        Inorm = np.reshape(Inorm.T, (h, w, 3)).astype(np.uint8)  
        
        # Create H Image
        H = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,0], axis=1).dot(np.expand_dims(C2[0,:], axis=0))))
        H[H>255] = 254
        H = np.reshape(H.T, (h, w, 3)).astype(np.uint8)
        #Create E Image
        E = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,1], axis=1).dot(np.expand_dims(C2[1,:], axis=0))))
        E[E>255] = 254
        E = np.reshape(E.T, (h, w, 3)).astype(np.uint8)
        
        return (Inorm, H, E,process)
        
    except np.linalg.LinAlgError:
        print("Warning: Eigenvalue decomposition failed. Skipping stain separation.")
        process = False
        return img, img, img,process # Return original image for all three outputs