function [ data ] = Fill5( data,fill_id)
% Fill5 is used to carry on linear interpolation
% 1 -99 -99 -99 -99 1can be filled by this function
[m,n]=size(data);
for j=1:n
    for i=2:m-5
        if(data(i,j)==fill_id&&data(i+1,j)~=fill_id)
            data(i,j)=(data(i-1,j)+data(i+1,j))/2;
        else if(data(i,j)==fill_id&&data(i+1,j)==fill_id&&data(i+2,j)~=fill_id)
                temp=(data(i-1,j)-data(i+2,j));
                data(i,j)=data(i-1,j)-temp/3;
                data(i+1,j)=data(i,j)-temp/3;
            else if(data(i,j)==fill_id&&data(i+1,j)==fill_id&&data(i+2,j)==fill_id...
                        &&data(i+3,j)~=fill_id)
                    temp=(data(i-1,j)-data(i+3,j));
                    data(i,j)=data(i-1,j)-temp/4;
                    data(i+1,j)=data(i,j)-temp/4;
                    data(i+2,j)=data(i+1,j)-temp/4;
                else if(data(i,j)==fill_id&&data(i+1,j)==fill_id&&data(i+2,j)==...
                            fill_id&&data(i+3,j)==fill_id&&data(i+4,j)~=fill_id)
                        temp=(data(i-1,j)-data(i+4,j));
                        data(i,j)=data(i-1,j)-temp/5;
                        data(i+1,j)=data(i,j)-temp/5;
                        data(i+2,j)=data(i+1,j)-temp/5;
                        data(i+3,j)=data(i+2,j)-temp/5;
                    else if((data(i,j)==fill_id&&data(i+1,j)==fill_id&&...
                                data(i+2,j)==fill_id&&data(i+3,j)==fill_id...
                                &&data(i+4,j)==fill_id&&data(i+5,j)~=fill_id))
                            temp=(data(i-1,j)-data(i+5,j));
                            data(i,j)=data(i-1,j)-temp/6;
                            data(i+1,j)=data(i,j)-temp/6;
                            data(i+2,j)=data(i+1,j)-temp/6;
                            data(i+3,j)=data(i+2,j)-temp/6;
                            data(i+4,j)=data(i+3,j)-temp/6;
                        end
                    end
                end
            end
        end
    end
end

                    
        

end

