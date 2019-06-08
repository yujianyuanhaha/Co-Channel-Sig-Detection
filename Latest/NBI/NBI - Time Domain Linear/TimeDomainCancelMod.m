function [y_lms, y_rls, y_lms_dd, y_rls_dd] = TimeDomainCancelMod(r, z, pulse, N, Nps, sps, bits, bitsPerSym, eBW, delta, lambda_inv, SequenceLength)
        MF = conv(r, pulse);
        
        init_idx = findLoc(bits, bitsPerSym, SequenceLength, sps, Nps, r, z);
        x = MF(init_idx:4:end);
        
        %x = MF(1:4:end);
        
        w_lms(:,1) = [zeros(1,Nps)]';
        w_rls(:,1) = [zeros(1,Nps)]';
        w_lms_dd(:,1) = [zeros(1,Nps)]';
        w_rls_dd(:,1) = [zeros(1,Nps)]';
        
        P = 10*eye(Nps);
        P_dd = 10*eye(Nps);
        
        lenX = length(x)
        if lenX <= N+8
            x = [x,zeros(1,N+8-lenX)];
        end
            
        
        for i=1:N
            z = x(i:i+7).';


            % first we test assuming an infinitely long training sequence
            % lms
            %[y_lms(:,i), w_lms(:,i+1)] = LMS(z,w_lms(:,i),delta, temp_bits(i));
            [y_lms(:,i), w_lms(:,i+1)] = LMS(z,w_lms(:,i),delta, 2*bits(i)-1);
            % rls
            %[y_rls(:,i), w_rls(:,i+1),P] = RLS(z, w_rls(:,i), P, lambda_inv, temp_bits(i));
            [y_rls(:,i), w_rls(:,i+1),P] = RLS(z, w_rls(:,i), P, lambda_inv, 2*bits(i)-1);
            
            
            % decision directed approach
            if i <= SequenceLength
                %TrainingLMS = temp_bits(i);
                %TrainingRLS = temp_bits(i);
                TrainingLMS = 2*bits(i)-1;
                TrainingRLS = 2*bits(i)-1;
            else
                TrainingLMS = sign(real(w_lms_dd(:,i)'*z));
                TrainingRLS = sign(real(w_rls_dd(:,i)'*z));
            end
            
            % lms
            [y_lms_dd(:,i), w_lms_dd(:,i+1)] = LMS(z,w_lms_dd(:,i),delta,TrainingRLS);
            % rls
            [y_rls_dd(:,i), w_rls_dd(:,i+1),P_dd] = RLS(z, w_rls_dd(:,i), P_dd, lambda_inv, TrainingRLS);
            
        end
end