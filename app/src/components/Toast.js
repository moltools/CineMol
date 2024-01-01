import { ToastContainer } from "react-toastify";

const Toast = () =>  {
    return (
        <ToastContainer
            position="bottom-right"
            autoClose={5000}
            hideProgressBar={true}
            newestOnTop={false}
            closeOnClick={false}
            rtl={false}
            pauseOnFocusLoss={false}
            draggable={false}
            pauseOnHover={true}
            icon={({ type }) => {
                if (type === "success") return "🎉";
                if (type === "error") return "🚨";
                else return "ℹ️";
            }}
        />
    );
};

export default Toast;