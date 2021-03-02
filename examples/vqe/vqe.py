import qsharp

from Microsoft.Quantum.Samples.Chemistry.SimpleVQE import GetEnergyHydrogenVQE

from scipy.optimize import minimize

def run_program(var_params, num_samples) -> float:
    # run parameterized quantum program for VQE algorithm
    theta1, theta2, theta3 = var_params
    return GetEnergyHydrogenVQE.simulate(theta1=theta1, theta2=theta2, theta3=theta3, nSamples=num_samples)
 
def VQE(initial_var_params, num_samples):
    """ Run VQE Optimization to find the optimal energy and the associated variational parameters """
 
    opt_result = minimize(run_program,
                          initial_var_params,
                          args=(num_samples,),
                          method="COBYLA",
                          tol=0.000001,
                          options={'disp': True, 'maxiter': 200,'rhobeg' : 0.05})
 
    return opt_result
 
if __name__ == "__main__":
    # Initial variational parameters
    var_params = [0.001, -0.001, 0.001]
 
    # Run VQE and print the results of the optimization process
    # A large number of samples is selected for higher accuracy
    opt_result = VQE(var_params, num_samples=100)
    print(opt_result)
 
    # Print difference with exact FCI value known for this bond length
    fci_value = -1.1372704220924401
    print("Difference with exact FCI value :: ", abs(opt_result.fun - fci_value))
