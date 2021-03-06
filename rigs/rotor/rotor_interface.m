classdef rotor_interface < rtc_interface
    %ROTOR_INTERFACE  Interface to the real-time control device that controls
    %the rotor rig. 
    
    % Written by Oliver Frolovs (oliver.frolovs@bristol.ac.uk) 2019 based on the code by
    % David A.W. Barton (david.barton@bristol.ac.uk) 2015
    
    properties
        RIG_STATUS_OK                         = 0;
        RIG_STATUS_UNKNOWN_ERROR              = 2^0;
        RIG_STATUS_SETUP_NOT_COMPLETED        = 2^1;
        RIG_STATUS_PID_NUMERIC_ERROR          = 2^2;
        RIG_STATUS_SPEED_SAFETY_LIMIT_REACHED = 2^3;
        RIG_STATUS_MOTOR_VOLTAGE_CLIP_AT_MIN  = 2^4;
        RIG_STATUS_MOTOR_VOLTAGE_CLIP_AT_MAX  = 2^5;
    end
    
    methods
        function obj = rotor_interface()
            %ROTOR_INTERFACE  Construct an ROTOR_INTERFACE object.
        end
        
        function reboot(obj)
            obj.par.bbb_reboot = 1;
        end
        
        function clear_status(obj)
            obj.par.reset_status_flags = 1;
        end
        
        function reset_pid_error_terms(obj)
            obj.par.reset_pid_error_terms = 1;
        end
        
        function v = speed_safety_limit_reached(obj)   
            v = bitand(obj.par.status_flags, obj.RIG_STATUS_SPEED_SAFETY_LIMIT_REACHED);
        end
        
        % TODO clear_status (to replace reset_status and provide granularity)
        
        function status(obj)
            flags = obj.par.status_flags;
            
            if flags == obj.RIG_STATUS_OK
                disp("  * Rig status OK")
            else
                if bitand(flags, obj.RIG_STATUS_UNKNOWN_ERROR)
                    warning("Rig error: unknown error");
                end

                if bitand(flags, obj.RIG_STATUS_SETUP_NOT_COMPLETED)
                    warning("Rig error: rig set up after reboot has not been finished");
                end

                if bitand(flags, obj.RIG_STATUS_PID_NUMERIC_ERROR)
                    warning("Rig error: PID controller numerical error");
                end
                
                if bitand(flags, obj.RIG_STATUS_SPEED_SAFETY_LIMIT_REACHED)
                    warning("Rig error: safety stop (%.1f rpm) was activated at %.1f rpm", ...
                        obj.par.speed_safety_limit_rpm, obj.par.safety_triggered_speed_rpm);
                end
                
                if bitand(flags, obj.RIG_STATUS_MOTOR_VOLTAGE_CLIP_AT_MIN)
                    warning("Motor voltage was clipped at minimum value of %.1f", ...
                        obj.par.motor_min_voltage);
                end
                
                if bitand(flags, obj.RIG_STATUS_MOTOR_VOLTAGE_CLIP_AT_MAX)
                    warning("Motor voltage was clipped at maximum value of %.1f", ...
                        obj.par.motor_max_voltage);
                end
            end
        end
    end
    
end

