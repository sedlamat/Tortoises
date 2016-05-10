import cv2
import numpy as np
import math
import os

<<<<<<< HEAD
def display_img(img, window_name='img'):
    cv2.namedWindow(window_name,cv2.WINDOW_NORMAL)
    cv2.imshow(window_name,img)
=======

def display_img(img):
    cv2.namedWindow('img',cv2.WINDOW_NORMAL)
    cv2.imshow('img',img)
>>>>>>> 5ddb5b0b2027d488862efcc89d570f84aa79411a
    cv2.waitKey(0)

def nofunction(x):
    pass

def scharr_masks(img):
    kernel = np.array([[-3,0,3], [-10,0,10], [-3,0,3]])
    convx = cv2.filter2D(img,cv2.CV_64F,kernel)
    convy = cv2.filter2D(img,cv2.CV_64F,np.transpose(kernel))
    convx = np.absolute(convx)
    convy = np.absolute(convy)
    print convx
    convx = convx.astype(np.uint8)
    convy = convy.astype(np.uint8)
    conv = cv2.addWeighted(convx,0.5,convy,0.5,0)
    maxconv = np.max(conv,axis=2)
    display_img(maxconv)
    return convx, convy

tortoise_image_name = 'Tg31400.pnm'
directory_with_tortoise_images = os.path.expanduser('~') + \
                                 '/Images/Tortoises/'
path_to_tortoise_image = directory_with_tortoise_images + \
                         tortoise_image_name

img = cv2.imread(path_to_tortoise_image,1)
resize_koef = 512.0/max(img.shape)
img = cv2.resize(img,(0,0),fx=resize_koef,fy=resize_koef)
img = cv2.GaussianBlur(img,(0,0),2)

min_of_BGR_values = np.min(img,axis=2)
color_only = img - np.dstack((min_of_BGR_values,
                              min_of_BGR_values,
                              min_of_BGR_values))

scharr_masks(color_only)
scharr_masks(img)
display_img(img)


display_img(without_intesity)
display_img(min_from_bgr)
display_img(img)
display_img(without_intesity)
display_img(min_from_bgr)
canny_b = cv2.Canny(without_intesity[:,:,0],0,0)
canny_g = cv2.Canny(without_intesity[:,:,1],0,0)
canny_r = cv2.Canny(without_intesity[:,:,2],0,0)
canny_i = cv2.Canny(min_from_bgr,0,0)
display_img(canny_b)
display_img(canny_g)
display_img(canny_r)
display_img(canny_i)
canny = cv2.bitwise_or(canny_b,canny_g)
canny = cv2.bitwise_or(canny,canny_r)
canny = cv2.bitwise_or(canny,canny_i)
display_img(canny)
display_img(img)
#print without_intesity

sobel_derivatives(without_intesity[:,:,0],3)

bCx, bCy = roberts_mask(without_intesity[:,:,0])
gCx, gCy = roberts_mask(without_intesity[:,:,1])
rCx, rCy = roberts_mask(without_intesity[:,:,2])
bC = cv2.addWeighted(bCx,0.5,bCy,0.5,0)
gC = cv2.addWeighted(gCx,0.5,gCy,0.5,0)
rC = cv2.addWeighted(rCx,0.5,rCy,0.5,0)

bgrC = np.dstack((bC,gC,rC))
display_img(bgrC)
C = np.max(bgrC,axis=2)
display_img(C)

val,C_t = cv2.threshold(C,20,255,cv2.THRESH_BINARY)
C_tI = cv2.bitwise_and(C_t,canny_i)
display_img(canny_i)
display_img(C_t)
display_img(C_tI)

Ix, Iy = roberts_mask(min_from_bgr)

I = cv2.addWeighted(Ix,0.5,Iy,0.5,0)
display_img(I)

CI = cv2.add(C,I)
display_img(CI)
display_img(img)

mCx = np.round(np.max(bgrCx,axis=2))
sCx = np.round(np.std(bgrCx,axis=2))
print np.max(sCx)
mCx = mCx*255.0/np.max(mCx)
sCx = sCx*255.0/np.max(sCx)
print 'sdf', np.max(sCx)
mCx = mCx.astype(np.uint8)
sCx = sCx.astype(np.uint8)
print mCx
print mCx.dtype
display_img(mCx)
display_img(sCx)

#conv = cv2.normalize(conv,alpha=0.0,beta=1.0)

#display_img(conv)

#mg = cv2.resize(img,(0,0),fx=res_koef,fy=res_koef)
#isplay_img(img[:,:])
#mg = cv2.GaussianBlur(img,(0,0),5)
#ap = cv2.Laplacian(img[:,:,0],ddepth=0,ksize=9)
#display_img(lap)
#lap = cv2.equalizeHist(lap)
#lap = cv2.normalize(lap,alpha=0,beta=255,norm_type=cv2.NORM_MINMAX)
#print np.max(lap)
#display_img(lap)

# hsv_tor = cv2.cvtColor(tor,cv2.COLOR_BGR2HSV)
# #rint hsv_tor
# display_img(tor)
# display_img(hsv_tor[:,:,0])
# exit

# print 255/blue_coef

# print blue_coef,green_coef,red_coef

img = np.zeros((512,1024,3),dtype='uint8')
window_name = 'track'
cv2.namedWindow(window_name,cv2.WINDOW_NORMAL)
cv2.createTrackbar('b',window_name,100,255,nofunction)
cv2.createTrackbar('g',window_name,0,255,nofunction)
cv2.createTrackbar('r',window_name,0,255,nofunction)
cv2.createTrackbar('d_b',window_name,100,255,nofunction)
cv2.createTrackbar('d_g',window_name,0,255,nofunction)
cv2.createTrackbar('d_r',window_name,0,255,nofunction)
while(True):
    b = cv2.getTrackbarPos('b',window_name)
    g = cv2.getTrackbarPos('g',window_name)
    r = cv2.getTrackbarPos('r',window_name)
    d_b = cv2.getTrackbarPos('d_b',window_name)
    d_g = cv2.getTrackbarPos('d_g',window_name)
    d_r = cv2.getTrackbarPos('d_r',window_name)
    img[:,:512] = (b,g,r)
    img[:,512:] = (b+d_b,g+d_g,r+d_r)
    print 2*d_b+b,2*d_g+g,2*d_r+r
    #print math.log(1.0+(d_b)*(b+d_b)), math.log(1.0+(d_g)*(g+d_g)), math.log(1.0+(d_r)*(r+d_r))
    cv2.imshow(window_name,img)
    k = cv2.waitKey(1)
    if k == 27:
        break

def canny_trackbars(img_gray):
    ''' 
    Shows edges detected with Canny, viewed with trackbars. 
    '''
    
    window_name = 'Canny - for exit press ESC'
    cv2.namedWindow(window_name,cv2.WINDOW_NORMAL)
    cv2.createTrackbar('lowThresh',window_name,50,100000,nofunction)
    cv2.createTrackbar('highThresh',window_name,150,100000,nofunction)
    cv2.createTrackbar('maskSize',window_name,3,7,nofunction)
    cv2.createTrackbar('isL2grad',window_name,0,1,nofunction)
    cv2.createTrackbar('blurSig*10',window_name,16,800,nofunction)
    
    while(True):
        thresh1 = cv2.getTrackbarPos('lowThresh',window_name)
        thresh2 = cv2.getTrackbarPos('highThresh',window_name)
        mask_size = cv2.getTrackbarPos('maskSize',window_name)
        is_L2grad = cv2.getTrackbarPos('isL2grad',window_name)
        sigma = cv2.getTrackbarPos('blurSig*10',window_name)
        
        if mask_size % 2 == 0:
            mask_size += 1
        if mask_size == 1:
            mask_size = 3

        sigma /= 10.0
        k_size = int(math.ceil(6*sigma))
        if k_size % 2 == 0:
            k_size += 1
            
        img_blured = cv2.GaussianBlur(src=img_gray,
                                      ksize=(k_size,k_size),
                                      sigmaX=sigma,
                                      dst=None,
                                      sigmaY=sigma,
                                      borderType=cv2.BORDER_DEFAULT)                              
        edges = cv2.Canny(image=img_blured,
                          threshold1=thresh1,
                          threshold2=thresh2,
                          apertureSize=mask_size,
                          L2gradient=is_L2grad)
        cv2.imshow(window_name,edges)
        k = cv2.waitKey(1) & 0xFF
        if k == 27:
            break
        
    cv2.destroyAllWindows()
    
def gradient_trackbars(img_gray):
    ''' 
    Shows image gradient, using Gaussian blur and Sobel masks, viewed with trackbars. 
    '''
    
    window_name = 'Gaussian blur - for exit press ESC'
    cv2.namedWindow(window_name,cv2.WINDOW_NORMAL)
    cv2.createTrackbar('blurSig*10',window_name,16,800,nofunction)
    cv2.createTrackbar('maskSize',window_name,3,7,nofunction)
    
    while(True):
        mask_size = cv2.getTrackbarPos('maskSize',window_name)
        sigma = cv2.getTrackbarPos('blurSig*10',window_name)
        
        if mask_size % 2 == 0:
            mask_size += 1
        if mask_size == 1:
            mask_size = 3

        sigma /= 10.0
        img_blured = img_gray
        if sigma >= 0.1:
            img_blured = cv2.GaussianBlur(img_gray,(0,0),sigma)
        kernely = np.array([1,-1],np.float32)/2
        dsty = cv2.filter2D(img_blured,cv2.CV_64F,kernely)
        dsty = np.absolute(dsty)
        kernelx = np.array([[1,-1]],np.float32)/2
        dstx = cv2.filter2D(img_blured,cv2.CV_64F,kernelx)
        dstx = np.absolute(dstx)
        dst = cv2.addWeighted(dstx,0.5,dsty,0.5,0)
        #sobelx = cv2.Sobel(img_blured,cv2.CV_64F,1,0,ksize=mask_size)
        # sobelxabs = np.absolute(sobelx)
        # sobely = cv2.Sobel(img_blured,cv2.CV_64F,0,1,ksize=mask_size)
        # sobelyabs = np.absolute(sobely)
        # sobel = cv2.addWeighted((sobelxabs,0.5,sobelyabs,0.5,0)

        cv2.imshow(window_name,dst)
        k = cv2.waitKey(1) & 0xFF
        if k == 27:
            break
        
    cv2.destroyAllWindows()

# if __name__ == "__main__":
#     main()
