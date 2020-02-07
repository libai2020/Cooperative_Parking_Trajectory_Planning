function WriteBoundaryValues2()
global vehicle_TPBV_
delete('BV');
fid = fopen('BV', 'w');
for ii = 1 : size(vehicle_TPBV_,2)
    fprintf(fid, '%g 1  %f\r\n', ii, vehicle_TPBV_{1,ii}.x0);
    fprintf(fid, '%g 2  %f\r\n', ii, vehicle_TPBV_{1,ii}.y0);
    fprintf(fid, '%g 3  %f\r\n', ii, vehicle_TPBV_{1,ii}.theta0);
    fprintf(fid, '%g 4  %f\r\n', ii, vehicle_TPBV_{1,ii}.xtf);
    fprintf(fid, '%g 5  %f\r\n', ii, vehicle_TPBV_{1,ii}.ytf);
    fprintf(fid, '%g 6  %f\r\n', ii, vehicle_TPBV_{1,ii}.thetatf);
end
fclose(fid);