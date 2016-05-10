
import os
import math
import numpy as np
import time
import cv2
import copy


def display_img(img,window_name='img'):
    cv2.namedWindow(window_name,cv2.WINDOW_NORMAL)
    cv2.imshow(window_name,img)
    cv2.waitKey(0)

def get_gradient_scharr(img):
    kernel = np.array([[-3,0,3], [-10,0,10], [-3,0,3]])
    dx = cv2.filter2D(img,cv2.CV_32F,kernel)
    dy = cv2.filter2D(img,cv2.CV_32F,np.transpose(kernel))
    return dx, dy

def get_gradient_magnitude(img):
    dx, dy = get_gradient_scharr(img)
    if len(dx.shape) == 3:
        dx = np.max(dx,axis=2)
        dy = np.max(dy,axis=2)
    return cv2.magnitude(dx,dy)

def get_gradient_orientation(img):
    dx, dy = get_gradient_scharr(img)
    if len(dx.shape) == 3:
        dx = np.max(dx,axis=2)
        dy = np.max(dy,axis=2)
    return np.arctan2(dx,dy) * 180 / np.pi

def get_edges_combo(img):
    img = np.copy(img)
    img = cv2.GaussianBlur(img,(0,0),1)
    min_of_BGR = np.min(img,axis=2)
    color_only = img - np.dstack((min_of_BGR,
                                  min_of_BGR,
                                  min_of_BGR))
    display_img(color_only)
    mag_color_only = get_gradient_magnitude(color_only)
    mag_min_of_BGR = get_gradient_magnitude(min_of_BGR)
    mag = mag_color_only*mag_min_of_BGR
    max_val = np.max(mag)
    koef = max_val/255.0
    edges_mag = cv2.threshold(mag,25000,255,cv2.THRESH_BINARY)[1]
    edges_mag = edges_mag.astype(np.uint8)
    edges = cv2.Canny(min_of_BGR,0,50)
    edges = cv2.bitwise_and(edges,edges_mag)
    return edges

def get_quantized_orientations(edges_orientations):
    orientations = np.copy(edges_orientations)
    range_idx = 0
    range_val = (1,45,90,135)
    ranges = ((23,67),(68,112),(113,157),(158,22))
    in_range = (orientations >= 158) + (orientations <= 22)
    orientations[in_range] = range_val[range_idx]
    range_idx += 1
    for lower_bound, upper_bound in ranges[:3]:
        in_range = (orientations >= lower_bound) * (orientations <= upper_bound)
        orientations[in_range] = range_val[range_idx]
        range_idx += 1
    return orientations

def get_r_table(template_orientations,reference_point):
    quantizations = (1,45,90,135)
    yr = reference_point[0]
    xr = reference_point[1]
    r_table = []
    for quant_val in quantizations:
        orient_positions = (template_orientations==quant_val)
        positions = [tuple(pos) for pos in np.argwhere(orient_positions)]
        positions_diff = [ ( yr-y, xr-x ) for y,x in positions]
        r_table.append(positions_diff)
    return r_table
                    
def get_accumulator(r_table,img_orientations):
    #display_img(img_orientations)
    range_val = (1,45,90,135)
    accumulator = np.zeros(img_orientations.shape)
    #print accumulator.dtype
    assert len(range_val)==len(r_table)
    for quant_idx, quant_val_diffs in enumerate(r_table):
        print quant_idx,
        #print len(quant_val_diffs)
        orient_positions = (img_orientations==range_val[quant_idx])
        positions = [tuple(pos) for pos in np.argwhere(orient_positions)]
        #print len(positions)
        #print len(positions)/100*0.05*4/60
        for y,x in positions:
            t0 = time.time()
            for y_diff, x_diff in quant_val_diffs:
                yr = y_diff + y
                xr = x_diff + x
                a = None
                b = None
                # where both yr,xr and y,x satisfy y = a*x + b
                if x_diff==0:
                    a = np.inf
                else:
                    a = y_diff*1.0/x_diff
                    b = (y*xr - yr*x)*1.0/x_diff

                if abs(a) < 1:
                    for xx in xrange(0,accumulator.shape[1],8):
                        if xx > x+10 or xx < x-10:
                            yy = a*xx + b
                            yy_f = int(math.floor(yy))
                            yy_c = int(math.ceil(yy))
                            if yy_f >= 0 and yy_f < accumulator.shape[0]:
                                accumulator[yy_f,xx] += 1
                            if yy_c >= 0 and yy_c < accumulator.shape[0]:
                                accumulator[yy_c,xx] += 1
                elif abs(a)==np.inf:
                    xx = x
                    for yy in xrange(0,accumulator.shape[0],8):
                        if yy > y+10 or yy < y-10:
                            accumulator[yy,xx] += 1
                else:       
                    for yy in xrange(0,accumulator.shape[0],8):
                        if yy > y+10 or yy < y-10:
                            xx = (yy-b)/a
                            #print xx,a,b,yy,x,y,x_diff,y_diff
                            xx_f = int(math.floor(xx))
                            xx_c = int(math.ceil(xx))
                            if xx_f >= 0 and xx_f < accumulator.shape[1]:
                                accumulator[yy,xx_f] += 1
                            if xx_c >= 0 and xx_c < accumulator.shape[1]:
                                accumulator[yy,xx_c] += 1
            #print((time.time()-t0))
            #print np.sum(img_orientations>0)
            #print img_orientations.size
        #print quant_idx
    return accumulator

def get_edges_orientations(img, img_edges):
    edges_orientations = get_gradient_orientation(img)
    edges_orientations += (edges_orientations<0)*180
    edges_orientations = edges_orientations.astype(np.uint8)
    edges_orientations = get_quantized_orientations(edges_orientations)
    edges_orientations = cv2.bitwise_and(edges_orientations,img_edges)
    return edges_orientations

def get_rotated_r_tables(r_table,rotations):
    rotated_r_tables = []
    for idx, angle in enumerate(rotations):
        print angle
        rotated_r_table = copy.deepcopy(r_table)
        for iii in xrange(idx-1):
            last_angle = rotated_r_table.pop()
            rotated_r_table.insert(0,last_angle)
        cs = math.cos(angle*np.pi/180)
        sn = math.sin(angle*np.pi/180)
        for idx0, quant_val_diffs in enumerate(rotated_r_table):
            for idx1, (y_diff, x_diff) in enumerate(quant_val_diffs):
                y = int(round(cs*y_diff + sn*x_diff))
                x = int(round(-sn*y_diff + cs*x_diff))
                rotated_r_table[idx0][idx1] = (y,x)
        rotated_r_tables.append(rotated_r_table)
    return rotated_r_tables
        
def main():
    
    directory = os.path.expanduser('~') + r'/Images/Generalized_Hough_Transform/'
    tortoise_images = os.listdir(directory)
    print tortoise_images
    return
    img_file_name = directory + 'tortoise3.jpg'
    ght_template_file_name = directory + 'ght_template.bmp'
    ght_template_edges_file_name = directory + 'ght_template_edges.bmp'
    img = cv2.imread(img_file_name)
    # M = cv2.getRotationMatrix2D((img.shape[1]/2,img.shape[0]/2),45,1)
    # img = cv2.warpAffine(img,M,(0,0),borderMode=cv2.BORDER_TRANSPARENT)

    ght_template = cv2.imread(ght_template_file_name)
    ght_template_edges = cv2.imread(ght_template_edges_file_name)
    resize_koef = 256.0/max(img.shape)
    img = cv2.resize(img,(0,0),fx=resize_koef,fy=resize_koef)
    img_edges = get_edges_combo(img)
    
    ght_template = ght_template[:,:,0]
    ght_template_edges = ght_template_edges[:,:,0]
    ref_point = tuple(np.argwhere(ght_template_edges == 127)[0])
    ght_template_edges = cv2.threshold(ght_template_edges,
                                       200,
                                       255,
                                       cv2.THRESH_BINARY)[1]
    template_edges_orientations = get_edges_orientations(ght_template,
                                                         ght_template_edges)
    img_edges_orientations = get_edges_orientations(img,img_edges)
    
    r_table = get_r_table(template_edges_orientations,
                          ref_point)

    rotations = [angle for angle in xrange(0,360,45)]
    print rotations
    rotated_r_tables = get_rotated_r_tables(r_table,rotations)
    accu_results = []
    for rotation_idx, angle in enumerate(rotations):
        print angle
        accumulator = get_accumulator(rotated_r_tables[rotation_idx],
                                      img_edges_orientations)
        accumulator = accumulator*255/np.max(accumulator)
        accumulator = accumulator.astype(np.uint8)
        kernel_size = 17
        kernel = np.ones((kernel_size,kernel_size))/(kernel_size*kernel_size)
        accumulator_mask = cv2.threshold(accumulator,1,1,cv2.THRESH_BINARY)[1]
        display_img(accumulator_mask)
        accumulator = cv2.filter2D(accumulator,cv2.CV_64F,kernel)
        accumulator_mask = cv2.filter2D(accumulator_mask,cv2.CV_64F,kernel)
        accumulator *= accumulator_mask
        min_val, max_val, min_coor, max_coor = cv2.minMaxLoc(accumulator)
        accu_results.append((max_val,angle,max_coor))
        print
    print accu_results
    accu_results = sorted(accu_results,reverse=True)
    print accu_results
    y_max = accu_results[0][2][1]
    x_max = accu_results[0][2][1]
    print 'angle', angle
    img[y_max-5:y_max+5,x_max-5:x_max+5] = (0,0,255)
    display_img(img)

if __name__ == '__main__':
    main()
