import ReactDOM from "react-dom";

import "./index.css";
import Header from "./components/Header";
import WorkSpace from "./components/WorkSpace";

export default function App() {
    return (
        <div className="app">
            <Header />
            <WorkSpace />
        </div>
    );
};

ReactDOM.render(<App />, document.getElementById("root"));