pro example

;create catalog with angle and flux constrain

COMPILE_OPT IDL2 


Nmax=file_lines('NVSSCatalog_v2.txt')


pi=3.14159265359

Nside=64
N64=12*64*64
outmap=fltarr(1,N64)
mapcount=fltarr(1,N64)





openr,lun,'NVSSCatalog_v2.txt',/get_lun


Numbercounts=0

read_fits_map, 'any_fits_map_just_need_the_header.fits',mask, hdr, xhdr, /silent

for i=0,Nmax/2-1 do begin


readf,lun,t1,t2,t3,t4,t5,t6,t7
readf,lun,p1,p2,p3

ra=ten(t1,t2,t3)
dec=ten(t4,t5,t6)

GLACTC,ten(t1,t2,t3),ten(t4,t5,t6) , 2000, gl, gb, 1 
phi=gl*pi/180
theta=pi/2-gb*pi/180

ang2pix_ring, Nside , theta, phi, ipring


if   (t7 gt 15)   then begin
outmap[ipring]++
Numbercounts++
endif
endfor


maskcount=0
pixelmean=0
for i=0,N64-1 do begin
pix2ang_ring, Nside,i, theta, phi
gb=180*(pi/2-theta)/pi
gl=phi*180/pi
GLACTC, ra, dec,  2000, gl, gb, 2



if  (mask[i] ne 0 )  then begin
pixelmean=pixelmean+outmap[i]
maskcount++
endif
end

print,'mask count',maskcount

pixelmean=pixelmean/maskcount

print,'mean',pixelmean


for h=0,N64-1 do begin
pix2ang_ring, Nside,h, theta, phi
gb=180*(pi/2-theta)/pi
gl=phi*180/pi
GLACTC, ra, dec,  2000, gl, gb, 2



outmap[h]=(outmap[h]-pixelmean)/pixelmean

end


write_fits_map, 'NVSS-V6-15mJy-N64.fits', outmap, Ordering='Ring'


close,/all



print,'*******'
print,'*Done!*'
print,'*******'

end
