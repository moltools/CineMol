import React, { useEffect } from "react";
import ReactDOM from "react-dom/client";
import { toast } from "react-toastify";

import "./index.css";
import Header from "./components/Header";
import Toast from "./components/Toast";
import WorkSpace from "./components/WorkSpace";

export default function App() {

    useEffect (() => {
        toast.info(
            "Welcome to CineMol!", 
            { autoClose: false } 
        );
    }, []);

    return (
        <div className="app">
            <Header />
            <WorkSpace />
            <Toast />
        </div>
    );
};

const root = ReactDOM.createRoot(document.getElementById("root"));

root.render(<App />);