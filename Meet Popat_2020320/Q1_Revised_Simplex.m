 %  Question 2 input
 % c = [-1,-1,-1,-1];
 % b = [1,1,1,1,1,1,1];
 % a = [1,1,1,0;0,1,0,1;0,0,1,1;1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];

 % Sample input
 a=[-1,1;1,1;2,5];
 b= [11,27,90];
 c = [4,6];


% Question 3 input 
  % c = [1,1,1,1];
  % b = [1,1,1,1,-1,-1,-1,-1];
  % a = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1;-1,-1,0,0;0,-1,-1,0;0,-1,0,-1;0,0,-1,-1];
  revised(a,b,c)

% Find only min
% All constrains st <=
function revised(a,b,c)
    % m x n matrix
    n=length(c);
    m=length(b);    

    j=max(abs(c));
    count=n;non_basic=1:n;basic=zeros(1,m);

    % for max comment out this
    % % c=-c;


       for i=1:m
            if b(i)<0
                a(i,:)=-a(i,:);
                b(i)=-b(i);
            end
            count=count+1;
            c(count)=0;
            a(i,count)=1;
            basic(i)=count;

        end
        A=[-c;a];
        B_inv=eye(m+1,m+1);     %identity matrix
        B_inv(1,2:m+1)=c(basic);   %inverse matrix
        x_b=B_inv*[0; b'];      %xb
        % iterations
        flag=0;
        count=0;
        curr=0;
        while(flag~=1)
            [s,t]=min(B_inv(1,:)*A(:,non_basic));
            y=B_inv*A(:,non_basic(t));
            count=count+1;
            if(any(y(2:m+1)>0))
                if count>1 && curr==x_b(1)
                    flag=1;
                else
                    curr=x_b(1);
                    if(s>=0)
                        flag=1;

                        for i=1:n
                            found=0;
                            for j=1:m
                               if basic(j)==i
                                    fprintf('x%u = %d\n',i,x_b(1+j)); %print basic variables
                                    found=1;
                                end
                            end
                        end
                    fprintf('The optimal value of the objective function=%d\n',x_b(1));
                    if found==0
                        fprintf('%u = %d\n',i,0);
                    end
                    else
                        u=10*j;
                        for i=2:m+1
                            if y(i)>0
                                if (x_b(i)/y(i))<u
                                    u=(x_b(i)/y(i));
                                    v=i-1;
                                end
                            end
                        end
                        temp=basic(v);
                        basic(v)=non_basic(t);
                        non_basic(t)=temp;      %Replacing values
                        E=eye(m+1,m+1);         %identity matrix
                        E(:,1+v)=-y/y(1+v);
                        E(1+v,1+v)=1/y(1+v);
                        B_inv=E*B_inv;
                        x_b=B_inv*[0; b'];
                    end
                end
            end
        end
end
